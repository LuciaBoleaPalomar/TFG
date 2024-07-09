
/* DESCRIPCIÓN DEL PROGRAMA:

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// --------------------- NÚMEROS ALEATORIOS | START -----------------------------
#define PI 3.14159265
#define NormaRanu (2.3283063671E-10F)
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;

float RapuanoRandom(void){ // Genera un número aleatorio entre 0 y 1
    float r;
    ig1=ind_ran-24; // Modificamos los índices cada vez que generamos un número
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=(irr[ig1]+irr[ig2]);
    ir1=(irr[ind_ran]^irr[ig3]); // Número random generado (entre 0 y 2^32-1)
    ind_ran++; // Cambiamos la posición base para el siguiente número random
    r=ir1*NormaRanu; // Divido entre NormaRanu para obtener un número entre [0,1)
    return r; 
}

/* Esta función sólo se ejecutará 1 única vez, al principio del programa, ya que
sólo sirve para inicializar el generador RapuanoRandom(). */ 
void IniRapu(void){ 
    double max;
    srand(time(NULL));
    max=4294967296;    //2^32
    for (int i=0; i<256; i++)
        irr[i]=max*(rand()/((double)RAND_MAX+1));

    ind_ran=ig1=ig2=ig3=0;
    for (int i=0; i<1024; i++)
        RapuanoRandom();
}

double gaussian(double x, double mean, double std_dev) {
    return exp(-0.5 * pow((x - mean) / std_dev, 2)) / (std_dev * sqrt(2 * PI));
}

// Generar un número aleatorio con distribución normal utilizando el método de Box-Muller
double rand_normal(double mean, double std_dev) {
    double u1 = RapuanoRandom();
    double u2 = RapuanoRandom();
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
    return z0 * std_dev + mean;
}
// --------------------- NÚMEROS ALEATORIOS | END -------------------------------



// PARÁMETROS GLOBALES:
int Ncells = 144; // number of cells
double tf = 1000; // tomamos como tiempo inicial t0 = 0
double h = 0.01; // diferencial de tiempo (el tiempo es discreto)

/* - FUNCIÓN A RESOLVER:
    f(t, y(t), y(t-tau))
    *xpunto: Array que devolverá las tasas de cambio de concentración en el tiempo t (evolución de x[time][cell]).
    i: iteración de RK2 (sirve para acceder a la solución en el tiempo desfasado tau)
    Ahora introduzco el retardo intercelular tau_c y el array x_tau con las concentraciones promedio de las células vecinas
*/
void f(double xpunto[Ncells], double x_tau[][Ncells], double x[Ncells], int i, double *tau, double beta, double K, double alpha, double Km, double V, double gam, double eps, double alph[][Ncells],double tau_c, double beta_c, double Kc){
    int i_delayed=0, i_couple=0;
    for(int j=0; j<Ncells; j++){
        i_delayed = i-round(tau[j]/h); // Índice correspondiente al tiempo retrasado de la célula j. Necesario para acceder a los distanciaes anteriores de x[][j].
        i_couple = i-round(tau_c/h); // Índice correspondiente al tiempo retrasado por el acoplo. Necesario para acceder a los distanciaes promedios anteriores de las células vecinas x[][neighbor].
        double x_ext=0;
        for (int k=0; k<Ncells; k++){
            // Media pesada con alph[j][k]
            x_ext+=alph[j][k]*x_tau[i_couple][k];
        }
        xpunto[j] = alpha/(1+pow(x_tau[i_delayed][j]/K,beta))*((1-eps)+eps/Kc/(1+pow(x_ext/Kc,beta_c)))-V*x[j]/(Km+x[j])-gam*x[j];
    }
}

/* - RK: Algoritmo Runge-Kutta de orden 2 para sistemas generales de N ecuaciones con N incógnitas.
    *f:    Puntero a la función vectorial
    t0:    Tiempo inicial
    tf:    Tiempo final
    x0:    Condición inicial
    N:     Dimensión del sistema (número de ecuaciones y de incógnitas); dimensión de los arrays x0, t, x, y
    x[time][cell]:  Array que contiene la evolución temporal de la función incógnita
*/
void RK2(double x[][Ncells], double *tau, double tau_max, double beta, double K, double alpha, double Km, double V, double gam, double eps, double alph[][Ncells], double tau_c, double beta_c, double Kc){
    int n_tau = tau_max/h; // Número de puntos almacenados con las condiciones iniciales

    double k1[Ncells];
    double k2[Ncells];
    double kAux[Ncells];

    /* Me coloco en el paso temporal correspondiente a t=tau porque se supone que hasta ahí ya tengo las 
    concentraciones guardadas en x. Quiero calcular el resto a partir de t=tau. */
    for (int i=n_tau; i<2*n_tau; i++){ // Usar el punto i-esimo para calcular el (i+1)-esimo.
        f(k1, x, x[i], i, tau, beta, K, alpha, Km, V, gam, eps, alph, tau_c, beta_c, Kc); // k1 = f(t_i,y_i)
        /* Le paso a la función:
            k1:   el array donde quiero que me guarde el cambio en las concentraciones (output)
            x:    matriz con todas las concentraciones de cada célula y en cada instante de tiempo
            x[i]: array de las concentraciones de las Ncells células para un cierto instante de tiempo correspondiente al índice i
            (el resto de argumentos son los que necesitamos en la fórmula) */
        // Con esta línea obtengo las tasas de cambio de las Ncells células en el instante i.

        /* Para cada célula, actualizo la concentración en el instante i en función de las tasas de cambio que 
        acabamos de guardar en k1.
        x[i][j] = concentración anterior
        k1[j]*k = (tasa de cambio en la concentración de la célula j) * dt      */
        for (int j=0; j<Ncells; j++){ 
            kAux[j] = x[i][j] + k1[j]*h;
        }
        f(k2, x, kAux, i+1, tau, beta, K, alpha, Km, V, gam, eps, alph, tau_c, beta_c, Kc); // k2 = f(t_i+h,y_i+h*k1)
        /* Le paso a la función:
            k2: el array donde quiero que me guarde este otro cambio en las concentraciones (output)
            x: matriz con todas las concentraciones de cada célula y en cada instante de tiempo
            kAux: array de las concentraciones en el tiempo correspondiente al índice i+h (línea 107) */
        // En este caso quiero calcular la tasa en el instante que corresponde al índice i+1.

        /* Calculo la concentración de todas las células en el instante i+1: */
        for (int j=0; j<Ncells; j++){
            x[i+1][j] = x[i][j] + h*0.5*(k1[j]+k2[j]);
        }
    }// bucle del instante i
}

int main(){
    IniRapu();
    
    // Parámetros de las células individuales
    double beta = 2.0;
    double K = 0.15;
    double alpha = 1.0;
    double Km = 50;
    double V = 1.0;
    double gam = 1.0; // Normalización implícita respecto al eje temporal
    
    // Asignación de los delays de cada célula
    double tau[Ncells];
    /*for (int j=0; j<Ncells; j++){
        tau[j]=3;
    }*/
    /*for (int j=0; j<Ncells; j++){
        if (j%2==0)
        tau[j]=3.3;
        else tau[j]=2.7;
    }*/
    /*for (int j=0; j<Ncells; j++){
        tau[j] = round((3.0+RapuanoRandom())*100.0)/100.0; // precisión máxima de cada delay: 2 decimales
        printf("\ntau[%d] = %f", j, tau[j]);
    }*/
    
    double l = 3;
    double media = 3;
    double desviacion_estandar = 0.05;
    const int steps = 10000; // Número de pasos de Metropolis
    for (int i=0; i<Ncells; i++) {
        for (int j = 0; j < steps; j++) {
            double l_new = l + rand_normal(0, desviacion_estandar); // Proponer un nuevo valor
            double acceptance_ratio = gaussian(l_new, media, desviacion_estandar) / gaussian(l, media, desviacion_estandar);
            if (acceptance_ratio >= 1.0 || RapuanoRandom() < acceptance_ratio) {
                l = l_new; // Aceptar el nuevo valor
            }
        }
        tau[i] = l; // Guardar el valor generado
    }
    
    
    // Parámetros de la interacción
    double eps = 1.0; // Interacción de la interacción (entre 0 y 1)
    double tau_c = 3; // Delay intercelular (siempre mayor que cualquier otro tau)
    double beta_c = 1.0;
    double Kc = 0.15;
    

    // --------------------- DEPENDENCIA CON LAS ARISTAS | START -------------------------------
    // Creo matrices para almacenar las distancias y los coeficientes de acoplo de cada célula
    double dist[Ncells][Ncells];
    double alph[Ncells][Ncells];
    for (int i=0; i<Ncells; i++) {
        for (int j=0; j<Ncells; j++) {
            dist[i][j] = -1; // Inicializo a -1 para ver mejor si hay errores
            alph[i][j] = 0;
        }
    }

    // Leer el fichero línea por línea
    FILE* fich = fopen("20240627.0.12x12.distances.txt", "r");
    if (fich == NULL) {
        printf("No se pudo abrir el archivo dist.\n");
    }
    char linea[100]; // Suponiendo que cada línea tiene como máximo 100 caracteres
    while (fgets(linea, sizeof(linea), fich)) {
        int ind_celula, ind_vecina;
        double distancia;
        sscanf(linea, "[%d, %d, %lf]", &ind_celula, &ind_vecina, &distancia);
        dist[ind_celula][ind_vecina] = distancia;
    }
    fclose(fich);

    // Plotteo las distancias leídas para ver si son correctas.
    for (int i=0; i<Ncells; i++) {
        printf("Celula %d:\n", i);
        for (int j=0; j<Ncells; j++) {
            if (dist[i][j] != -1) {
                printf("Vecina: %d \t Distancia: %.15lf\n", j, dist[i][j]);
            }
        }
        printf("\n");
    }
    
    // Calculo los coeficientes de acoplo, dependientes de las distancias
    for (int i=0; i<Ncells; i++){
        double suma_distancias = 0;
        for (int j=0; j<Ncells; j++){
            if (dist[i][j]!=-1){
                suma_distancias+=dist[i][j];
            }
        }
        for (int j=0; j<Ncells; j++){
            if (dist[i][j]!=-1){
                alph[i][j] = dist[i][j]/suma_distancias;
            }
        }
        /*double compruebo = 0;
        for (int j=0; j<Ncells; j++){
            compruebo+=alph[i][j];
        }
        printf("%d, %lf", i, compruebo); */      //todo ok
    }
    // --------------------- DEPENDENCIA CON LAS ARISTAS | END -------------------------------
    


    // ----------------------- GUARDO INFO DEL TEJIDO | START -----------------------
    // Fichero con la evolución temporal de las CONCENTRACIONES 
    FILE *fich_x;
    fich_x=fopen("20240709.0.12x12.evol_temp_sigma0.05_tauc3-repe.txt","w");
    fprintf(fich_x,"#VARIABLES DEL SISTEMA:\n");
    fprintf(fich_x,"#beta = %.3f\t alpha = %.3f\t K = %.3f\t V = %.3f\t Km = %.3f\t gamma = %.3f\t\n", beta, alpha, K, V, Km, gam);
    fprintf(fich_x,"#t\t");
    for(int j=0; j<Ncells; j++){
        fprintf(fich_x,"x[%d]\t",j); // En cada columna tengo una célula
    }
    fprintf(fich_x,"\n");
    
    /*Fichero con la evolución temporal del PARÁMETRO DE ORDEN ------------------------------------------------
    FILE *fich_R;
    fich_R=fopen("20240706.1.5.12x12.order_R_sigma0_tauc3.txt","w");
    fprintf(fich_R,"#VARIABLES DEL SISTEMA:\n");
    fprintf(fich_R,"#beta = %.3f\t alpha = %.3f\t K = %.3f\t V = %.3f\t Km = %.3f\t gamma = %.3f\t\n", beta, alpha, K, V, Km, gam);
    //fprintf(fich_R,"#t_order\tR\n");
    fprintf(fich_R,"#tau_c\tR\n");
    */ // ----------------------- GUARDO INFO DEL TEJIDO | END -----------------------

    //Búsqueda del tau máximo
    double tau_max = tau_c;
    for (int i=0; i<Ncells; i++){
        if (tau[i]>tau_max)
            tau_max = tau[i];
    }
        
    // Variables para integración
    int n_tau = tau_max/h; // Número de puntos almacenados con las condiciones iniciales
    int iter_max = ceil(tf/tau_max); // Número de iteraciones de RK2 (por intervalos de tiempo de longitud tau_max) necesarias para resolver el sistema
        
    double t[2*n_tau+1]; // Almacenamos solo el intervalo actual y el anterior
    double x_cells[2*n_tau+1][Ncells];
        
    //INICIALIZACIÓN DEL SISTEMA: entre -tau y 0, cada célula toma distancia constante x0 asignado aleatoriamente
    double aux=0;
    for (int j=0; j<Ncells; j++){
        aux = RapuanoRandom(); // Condiciones iniciales aleatorias entre 0 y 1
        for (int i=0; i<=n_tau; i++){
            x_cells[i][j]=aux;
            /*x_cells[i][0]=0.1;
            x_cells[i][1]=0.7;*/
        }
    }


    // Pasar las concentraciones de las condiciones iniciales al fichero
    for(int i=0; i<=n_tau; i++){ // Escribo los n_tau pasos temporales (de momento sólo conozco estas concentr.)
        t[i]=i*h-tau_c; // Para mostrar los tiempos desde -tau_c hasta t=0
        fprintf(fich_x,"%.3f\t",t[i]); // Escribo los tiempos "antes del 0" (con el 0 incluido). En cada fila tengo un dt
        for(int j=0; j<Ncells; j++){
            fprintf(fich_x,"%.3f\t",x_cells[i][j]); // Escribo las concentraciones de las Ncells para cada dt
        }
        fprintf(fich_x,"\n");
    } // ------------------------------------------------------------------------------------------------------
        
    /* Hasta ahora tengo sólo la mitad (+ el elemento correspondiente al 0) de los elementos de x_cells y t 
    (hasta el índice n_tau, siendo el total de elementos 2*n_tau+1). Sólo he dado distanciaes a t desde -tau a 0, 
    y las concentraciones correspondientes son las concentraciones iniciales (distancia constante desde -tau a 0). */


        
        
    // Defino VARIABLES para calcular el PARÁMETRO DE ORDEN
    int iter_term = ceil(0.1*tf/tau_c); // TERMALIZACIÓN: iteración del total de simulacion a partir de la cual calculamos R
    double num=0, denom;
    double M, sum_M=0, sum_M2=0; // En sum_M almacenamos la suma de M para el promedio temporal. Íden con sum_M2 y la suma de cuadrados.
    double sum_xi[Ncells], sum_xi2[Ncells]; // En sum_xi almacenamos la suma de x para el promedio temporal. Íden con sum_xi2 y la suma de cuadrados.
    for(int j=0; j<Ncells; j++){
        sum_xi[j] = 0;
        sum_xi2[j] = 0;
    }
            
        
    //------------------- CÁLCULO DE LA DINÁMICA | START --------------------------------
        
    //INTEGRACIÓN DE LA ED POR INTERVALOS DE TIEMPO DE LONGITUD tau
    int t_steps=0; // Número de puntos usados para calcular el promedio temporal (contador de iteraciones global)
    int n_cycles=0; // Contador para calcular la frecuencia de oscilación

    for (int cont_iter=0; cont_iter<iter_max; cont_iter++){
        for (int i=0; i<=2*n_tau; i++){ // Actualizo el array de tiempo en la iteración cont_iter
            t[i]=cont_iter*tau_c+(i*h-tau_c);
            /* Los distanciaes de t van desde [cont_iter*tau]-tau a [cont_iter*tau]+tau (-tau a tau + el factor [ ])
                cont_iter = 0 ->  t[0]=-tau | t[2*n_tau]=tau
                cont_iter = 1 ->  t[0]=0    | t[2*n_tau]=2tau 
                cont_iter = 2 ->  t[0]=tau  | t[2*n_tau]=3tau    etc.    
            El array t abarca 2tau distanciaes, y con cada cont_iter, el inicio de t se desplaza tau instantes.    */ 

            /* El distancia en t[n_tau] para cont_iter=m, se repite en t[0] para cont_iter=m+1
            Por ejemplo, t[i] = 0 aparece en dos momentos diferentes:
            - Cuando cont_iter=0, en t[300]=t[n_tau]
            - Cuando cont_iter=1, en t[0] 
            */

            /* En la primera iteración de cont_iter en realidad reescribo los primeros n_tau elementos porque
            cuando he escrito las n_tau condiciones iniciales desde -tau_c hasta 0, ya he definido esas compo-
            nentes de t, pero da igual. Es solo un apunte.
            */
        }
            
        // Ejecutar RK2 para cada tramo para calcular x_cells en [cont_iter*tau, (cont_iter+1)*tau]
        RK2(x_cells, tau, tau_max, beta, K, alpha, Km, V, gam, eps, alph, tau_c, beta_c, Kc);

        // Escribo las concentraciones en el fichero -------------------------------------------------------------
        /* Hasta ahora había escrito las concentraciones en todas las células para los tiempos entre [-tau_c, 0].
        Ahora escribo las nuevas concentraciones que he guardado en la segunda mitad de x_cells. */
        for(int i=n_tau+1; i<=2*n_tau; i++){
            /*  Como ya he escrito desde -tau=t[0] a 0=t[n_tau] en las condiciones iniciales, ahora empiezo a escribir
            desde t[n_tau+1], hasta el último elemento calculado por RK2: t[2*n_tau]. 
            El elemento en t[2*n_tau] será t[n_tau] en el siguiente cont_iter, pero no se reescribirá porque empiezo
            a escribir desde t[n_tau+1].  */

            fprintf(fich_x,"%.3f\t",t[i]);
            for(int j=0; j<Ncells; j++){
                int index = round(i-tau[j]/h); // DUDA: No hace nada con esto?????
                fprintf(fich_x,"%.3f\t",x_cells[i][j]);
            }
            fprintf(fich_x,"\n");
        } // -----------------------------------------------------------------------------------------------------

        /* Sólo calculamos las cosas a partir de cierto número de iteraciones (termalización). 
            Concretamente a partir de cont_iter=34.                                                     */ 
        if(cont_iter>=iter_term){
            for (int i=n_tau+1; i<=2*n_tau; i++){ // Tomo los elementos de x_cells que acabo de calcular con RK2
                //printf("t_steps: %d\t", t_steps);
                /* En la próxima iteración de este for tendré otro índice i y querré calcular los 
                promedios añadiendo otro tiempo t[i], así que pongo a 0 las variables en las 
                que inserto la suma.   */
                M=0;
                denom=0;
                /* Sin embargo, sum_M, sum_M2, sum_xi y sum_xi2 no las pongo a 0 porque lo que
                hago es añadir a la suma que ya tenía los distanciaes de este nuevo t[i], para 
                luego dividir entre t_steps, que habrá aumentado una unidad.  */

                t_steps++; /* Esto tiene que aumentar en cada iteración del bucle de las i, porque 
                                para cada instante de tiempo (una vez he termalizado) se calcula num y 
                                denom, luego al hacer la media tengo que tener en cuenta todas las veces
                                que he sumado.   */ 

                // NUMERADOR
                for(int j=0; j<Ncells; j++){
                    M += x_cells[i][j]; // Cálculo de media de concentraciones en tiempo t[i]
                }
                M = M/Ncells;
                sum_M += M;
                sum_M2 += pow(M,2);
                num = sum_M2/t_steps-pow(sum_M/t_steps,2);
                    
                // DENOMINADOR
                for(int j=0; j<Ncells; j++){
                    sum_xi[j] += x_cells[i][j];
                    sum_xi2[j] += pow(x_cells[i][j],2);
                    denom += sum_xi2[j]/t_steps - pow(sum_xi[j]/t_steps,2);
                }
                denom = denom/Ncells;
                    
                // Escribo el parámetro de orden en el fichero -----------------------------------------
                // fprintf(fich_R,"%.3f\t%.3f\n", t[i], num/denom);
                //printf("t = %.3f, R = %.3f\n", t[i], num/denom);
                //printf("t[i]=%.2f\tR=%.3f\n", t[i], num/denom);
                // ------------------------------------------------------------------------------------

            } // Se cierra el for para este instante de tiempo i
        } // Se cierra el if(cont_inter>=iter_term).
            
        // Actualizar array de concentraciones
        for(int i=0; i<=n_tau; i++){
            for(int j=0; j<Ncells; j++){
                x_cells[i][j] = x_cells[i+n_tau][j];
            }
        }
    }  // Bucle for cont_iter

    //------------- CÁLCULO DE LA SOLUCIÓN NUMÉRICA ------------- END
    

    printf("All good.\n");
    fclose(fich_x);
    // fclose(fich_R);
    return 0;
}


