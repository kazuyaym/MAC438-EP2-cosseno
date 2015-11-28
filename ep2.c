/**************************************/
/*                                    */
/*   EP 2 - Programação Concorrente   */
/*                                    */
/*   Marcos Kazuya          7577622   */
/*   ep2.c                            */
/*                                    */
/**************************************/

#include "gmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#include <string.h>

#define MAX_SEMAFORO 1024
#define BITS_PRECISAO 1048576

typedef struct thread {
  pthread_t thread_id; 
  unsigned int id; /* ID da thread */
  unsigned int n; /* Sera a iteracao i da formula de Taylor da funcao cos(x) */
} THREAD;

/* String que sera usada para mudar o numero de caracteres decimais para serem impressas */
char buf[128], buff_parcial[128]; 

/* Argumentos recebidos na entrada */
int arg_qntThreads;    /* Quantidade de threads que serao criadas                   */
int arg_diffMenorQue;  /* Segundo argumento == f                                    */
int arg_valorMenorQue; /* Segundo argumento == m                                    */
int arg_precisao;      /* valor do terceiro argumento, que sera utilizado
                          de acordo com o valor do segundo argumento                */

int arg_debug;         /* argumento 5, d => 1, cc => 0                              */
int arg_sequencial;    /* argumento 5, s => 1, cc =>0                               */
mpf_t arg_cosx;   /* Valor de x para se calcular o cosseno em radiano               */
mpf_t x2;         /* arg_cosx ao quadrado para simplificar contar                   */

int continua;       /* Variavel que sera verificada na barreira
                       se o programa deve continuar ou parar                        */
int barreira;       /* Quantidade de vezes que as threads
                       passaram pela barreira                                       */

mpf_t termo_de_parada;  /* Numero que sera comparado para o programa terminar       */
mpf_t diferenca_termos; /* Caso o segundo argumento seja f, essa variavel
                           sera usada para fazer a comparacao se o
                           programa deve terminar                                   */
mpf_t cosx;         /* Valor final do cosseno(x) a cada iteracao                    */
mpf_t cosx_antigo;  /* Valor do cosseno da rodada anterior                          */
mpf_t termo;        /* Valor absoluto do termo de taylor calculado em cada iteracao */

sem_t mutex[MAX_SEMAFORO]; /* Serao usados arg_qntThreads semaforos                 */

/* Variaveis para a barreira */
int arrive;
pthread_mutex_t mutexBarreira;
pthread_cond_t condicaoBarreira;

/* Funcao em que cada thread executara                    */
/* Sera passado a struct THREAD como parametro de entrada */
void *cos_taylor(void *arg){
  THREAD *tinfo = arg; 
  unsigned int n2;

  /* Determina qual o proximo semaforo a ser liberado */
  int proxThread;

  mpf_t term1;
  mpf_init2(term1, BITS_PRECISAO);

  
   
  if (tinfo->id == arg_qntThreads - 1) proxThread = 0;
  else proxThread = tinfo->id + 1;

  while(continua) {
    if(tinfo->n) { /* Cada thread tera que calcular um termo x2/2n*(2n-1) */
      n2 = 2*tinfo->n;
      mpf_div_ui(term1, x2, n2*(n2-1));
    }
    else { /* Caso seja o primeiro, onde tinfo->n == 0, 
              temos isso, pois nao podemos dividir por zero */
      mpf_set_ui(term1, 1);
    }

    /* SECAO CRITICA */
    sem_wait(&mutex[tinfo->id]);

      /* calular o termo de taylor dessa iteracao */
      mpf_mul(termo, termo, term1);

      /* incrementar na variavel do cosseno */
      if(tinfo->n % 2) mpf_sub(cosx, cosx, termo);
      else mpf_add(cosx, termo, cosx);

    /* SECAO NAO CRITICA */
    sem_post(&mutex[proxThread]);

    /* Barreira, espera todos os processos a 
       chegarem neste estagio antes de proseguirem */
    pthread_mutex_lock(&mutexBarreira);
    arrive--; 
    if(arg_debug) printf("Thread %d chegou na barreira\n", tinfo->id);

    if ( arrive == 0 ) { 

      /* verificar se o termo achado nao 'e o de parada */
      if(arg_valorMenorQue && mpf_cmp(termo, termo_de_parada) < 0)
          continua = 0;
      else if(arg_diffMenorQue) {
        mpf_sub (diferenca_termos, cosx, cosx_antigo); /* subtrai o valor do cosseno novo, com o antigo */
        mpf_abs (diferenca_termos, diferenca_termos); /* Pega o valor absoluto da operacao anterior */
        if(mpf_cmp(diferenca_termos, termo_de_parada) < 0) continua = 0;
        else mpf_set (cosx_antigo, cosx);
      }

      if(arg_debug) gmp_printf(buff_parcial, arg_cosx ,cosx);
      
      barreira++;
      if(continua && arg_debug) printf("Rodada %d\n", barreira + 1);

      arrive = arg_qntThreads;
      pthread_cond_broadcast(&condicaoBarreira);
    }
    else
      pthread_cond_wait(&condicaoBarreira, &mutexBarreira); 
    pthread_mutex_unlock(&mutexBarreira);

    tinfo->n += arg_qntThreads; /* Proxima iteracao da somatoria na formula de Taylor */
  }
  mpf_clear(term1);
  return NULL;
}

void cos_taylor_sequencial() {
  unsigned int n2;
  int n = 1;

  mpf_t term1;
  mpf_init2(term1, BITS_PRECISAO);

  mpf_add(cosx, termo, cosx);
  gmp_printf(buff_parcial, arg_cosx, cosx);
  printf("Rodada %d\n", n+1);

  while(1) {
    n2 = 2*n;
    mpf_div_ui(term1, x2, n2*(n2-1));

    /* calular o termo de taylor dessa iteracao */
    mpf_mul(termo, termo, term1);

    /* incrementar na variavel do cosseno */
    if(n % 2) mpf_sub(cosx, cosx, termo);
    else mpf_add(cosx, termo, cosx);
    
    /* verificar se o termo achado nao 'e o de parada */
    if(mpf_cmp(termo, termo_de_parada) < 0) break;

    gmp_printf(buff_parcial, arg_cosx ,cosx);
    printf("Rodada %d\n", ++n+1);
  }

  mpf_clear(term1);
  return;
}

int main(int argc, char * argv[]){
  int i; /* Iterador usado na inicializacao do semaforo, e das threads */
  THREAD *thread_info;

  if (argc <= 4){
    printf ("Usage: %s <number of threads> <m|f> <precisao> <cos (x)> <opcional: d|s>\n", argv[0]);
    return EXIT_FAILURE;
  }
  
  arg_qntThreads = atoi(argv[1]);
  if(arg_qntThreads == 0) arg_qntThreads = sysconf(_SC_NPROCESSORS_ONLN);

  arg_valorMenorQue = 0;
  arg_diffMenorQue  = 0;
  if ( strcmp(argv[2], "m") == 0 )      arg_valorMenorQue = 1;
  else if ( strcmp(argv[2], "f") == 0 ) arg_diffMenorQue  = 1;

  arg_precisao = atoi(argv[3]);
  if( arg_precisao < 0) arg_precisao = arg_precisao * (-1);

  mpf_init2(arg_cosx, BITS_PRECISAO);
  mpf_set_d(arg_cosx, atof(argv[4]));

  arg_debug = 0;
  arg_sequencial = 0;
  if ( argc > 5){
    if ( strcmp(argv[5], "d") == 0 ) arg_debug = 1;
    if ( strcmp(argv[5], "s") == 0 ) arg_sequencial = 1;
  }
    
  /* inicializa a variavel onde sera guardada o valor do cosseno de x */
  mpf_init2(cosx, BITS_PRECISAO); 
  mpf_init2(termo, BITS_PRECISAO);
  mpf_init2(termo_de_parada, BITS_PRECISAO);
  if(arg_diffMenorQue) {
    mpf_init2(cosx_antigo, BITS_PRECISAO);
    mpf_init2(diferenca_termos, BITS_PRECISAO);
  }

  mpf_set_ui(termo_de_parada, 10); /* seta termo de parada em 10       */ 
  if(arg_precisao > 0) { /* Se caso arg_precisao for negativo teriamos */
                         /* 10^(-arg_precisao), que seria >= 10,       */
                         /* e o programa pararia na primeira rodada    */
    mpf_pow_ui(termo_de_parada, termo_de_parada, arg_precisao); /* Eleva termo de parada com arg_precisao. Ex 10^1000 */
    mpf_ui_div(termo_de_parada,1,termo_de_parada); /* Como queriamos 10^(-1000), fazemos 1/termo_de_parada */
  }

  /* Cria a String que sera passada como primeiro argumento no gmp_printf */
  arg_precisao *= 2;
  if(arg_precisao > 100000) arg_precisao = 100000; /* o maximo de caracteres que serao impressos sera 100k */
  if(arg_precisao < 10    ) arg_precisao = 10; /* o minimo de caracteres que vamos imprimir eh 10 */
  sprintf(buf, "Cosseno(%%.5Ff) = %%.%dFf\n", arg_precisao); /* puts string into buffer */
  sprintf(buff_parcial, "Valor parcial do Cosseno(%%.5Ff) = %%.%dFf\n\n", arg_precisao); /* puts string into buffer */
  
  mpf_set_ui(cosx, 0);
  mpf_set_ui(termo, 1);

  mpf_init2(x2, BITS_PRECISAO);
  mpf_mul(x2, arg_cosx, arg_cosx); /* x ao quadrado, para diminuir as contas nas threads */
  
  if(arg_debug || arg_sequencial) printf("Rodada 1\n");

  if(!arg_sequencial) {

    continua = 1;
    barreira = 0;

    /* Inicializa array da struct onde tem as informacoes de cada thread */
    thread_info = calloc(arg_qntThreads, sizeof(THREAD));
    if (!thread_info) {
      printf("Não foi possível alocar memória!\n");
      return EXIT_FAILURE;
    }

    /* Inicializa o semáforo */
    sem_init(&mutex[0],0,1); /* mutex[0] esta liberado */
    for(i = 1; i < arg_qntThreads; i++){ /* outros semaforos serao liberados apenas
                                quando a thread de id-1 sair da seccao
                                critica */
      if (sem_init(&mutex[i],0,0)) {
          printf("Erro ao criar o semáforo!\n");
          return EXIT_FAILURE;
      }
    }

    /* Inicializando mutex, condicao da barreira e a variavel arrive 
       que contem a qnt de threads que a barreira deve esperar */
    pthread_mutex_init(&mutexBarreira, NULL);
    pthread_cond_init(&condicaoBarreira, NULL);
    arrive = arg_qntThreads;

    /* cria cada thread, atribuindo 
       thread_info[i].id = i, que sera inalterado durante todo o processo
       thread_info[i].n = i, que sera alteraco a cada iteracao para saber qual
                             valor do termo de taylor sera calculado           */
    for(i = 0; i < arg_qntThreads; i++) {
      thread_info[i].id = i;
      thread_info[i].n = i;

      if ( pthread_create(&thread_info[i].thread_id, NULL, &cos_taylor, &thread_info[i]) ) {
        printf("Não foi possível criar a thread %d \n", thread_info[i].id);
        return EXIT_FAILURE;
      }

    }

    /* Espera todos os processos chegarem a esse ponto */
    for (i = 0; i < arg_qntThreads; ++i) {
      if ( pthread_join(thread_info[i].thread_id, NULL) ) {
        printf("Não foi possível dar join na thread %d \n", thread_info[i].id);
        return EXIT_FAILURE;
      }
    }
    printf("As threads passaram pela barreira %d vezes\n", barreira);
    free(thread_info);
  }
  else /* Calcula o valor de cosseno em sequencial */
    cos_taylor_sequencial();

  /*gmp_printf("Cosseno(%d) = %.100000Ff\n", arg_cosx ,cosx);*/
  gmp_printf(buf, arg_cosx ,cosx);

  mpf_clear(cosx);
  mpf_clear(termo);
  mpf_clear(termo_de_parada);
  mpf_clear(arg_cosx);
  mpf_clear(x2);
  if(arg_diffMenorQue) {
    mpf_clear(cosx_antigo);
    mpf_clear(diferenca_termos);
  }
  return EXIT_SUCCESS;
}