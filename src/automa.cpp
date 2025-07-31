#include "../include/automa.h"
#include <cmath>
#include <stdio.h>

Locazione::Locazione(TipoStato s){
  this->stato=s;
  this->f_ODE=nullptr;
  this->condizione=nullptr;
}

Automa::Automa(set<TipoStato>& stati){
  this->locazioni=set<Locazione>();
  for(unsigned long s : stati){
    this->locazioni.insert((const unsigned long)s);
  }
}

Automa::Automa(set<TipoStato> stati){
  this->locazioni=set<Locazione>();
  for(unsigned long s : stati){
    this->locazioni.insert((const unsigned long)s);
  }
}

bool Automa::AggiungiRegola(TipoStato stato, TipoInput input, TipoStato statoSucc, TipoOutput output){
  bool ret=false;
  if(this->locazioni.contains(stato) and this->locazioni.contains(statoSucc)){  //Se sono presenti gli stati nell'automa
    pair<TipoStato,TipoInput> chiave{stato,input};
    pair<TipoStato,TipoOutput> valore{statoSucc,output};
    this->transizioni.insert({chiave,valore});
    this->guardie[chiave]=pair<Condizione,FunzioneReset>({nullptr,nullptr});
    ret=true;
  }
  return ret;
}

pair<queue<TipoStato>,queue<TipoOutput>> Automa::ValutaInput(queue<TipoInput> stringaInput,TipoStato s0){
  queue<TipoStato> evoluzione;
  queue<TipoOutput> uscita;
  if(this->locazioni.contains(s0)){
    evoluzione.push(s0);
    TipoStato statoAttuale=s0;
    for(;not stringaInput.empty();stringaInput.pop()){
      TipoInput i=stringaInput.front();
      pair<TipoStato,TipoInput> coppia({statoAttuale,i});
      if(this->transizioni.contains(coppia)){
        pair<TipoStato,TipoOutput> valore=this->transizioni[coppia];
        evoluzione.push(valore.first);
        uscita.push(valore.second);
        statoAttuale=valore.first;
      }else break;
    }
  }
  return pair<queue<TipoStato>,queue<TipoOutput>>({evoluzione, uscita});
}

void Automa::ImpostaODE(TipoStato s, ODE f){
  if(this->locazioni.contains(s)){
    auto NodoLocazione=this->locazioni.extract(s);
    NodoLocazione.value().f_ODE=f;
    this->locazioni.insert(move(NodoLocazione));
  }
}

void Automa::ImpostaCondizioneLocazione(TipoStato s, Condizione f){
  if(this->locazioni.contains(s)){
    auto NodoLocazione=this->locazioni.extract(s);
    NodoLocazione.value().condizione=f;
    this->locazioni.insert(move(NodoLocazione));
  }
}

void Automa::ImpostaCondizioneGuardia(TipoStato s, TipoInput i, Condizione g, FunzioneReset r){
  pair<TipoStato,TipoInput> coppia({s,i});
  if(this->transizioni.contains(coppia)){
    this->guardie[coppia]=pair<Condizione,FunzioneReset>({g,r});
  }
}

pair<gsl_matrix*,queue<pair<double,TipoStato>>> Automa::Simulazione(double t0, gsl_vector* x0, TipoStato s0, double T, double h, MetodoODE mODE, MetodoInnesco mInnesco, unsigned passi){
  size_t n=x0->size, p=passi-1, NumeroCampioni=(size_t)floor(pow(10.0,log10(T)-log10(h)))+1;
  queue<pair<double,TipoStato>> evoluzioneAutoma=queue<pair<double,TipoStato>>();
  gsl_matrix* evoluzioneODE=gsl_matrix_calloc(n,NumeroCampioni);
  gsl_matrix_view innesco;
  struct InfoBaseSimulazione infoSimulazione;
  double istanteCondizione=0.0,t=t0;
  size_t indiceCondizione=0;
  bool innescoPronto=false,bloccante=true;
  TipoStato statoAttualeAutoma=s0;
  gsl_vector* statoInizialeODE=gsl_vector_calloc(n);
  gsl_vector_memcpy(statoInizialeODE,x0);
  
  size_t o=0;
  for(size_t indiceIstante=0; indiceIstante < NumeroCampioni-1 and o < 20;){
    ++o;
    innescoPronto=false;
    bloccante=false;
    auto locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
    
    printf("Campioni rimanenti: %lu\tt: %.10lf\n",NumeroCampioni-indiceIstante,t);
    size_t campioniRimanenti=NumeroCampioni-indiceIstante;
    
    //Preparazione innesco
    while(not innescoPronto){
      if(p > 0){
        //Calcolo innesco
        //gsl_matrix* calcoloInnesco=Innesco()
      }else{
        printf("p=0\n");
        bool verificaCondizione= locazioneAttuale->condizione ? locazioneAttuale->condizione(t,statoInizialeODE) : false;
        if(verificaCondizione){
          bloccante=true;
          for(auto transizione : this->transizioni){ //Se il ciclo for termina allora non ci sono transizioni abilitate per uscire dalla locazione, il sistema è bloccante
            if(transizione.first.first == statoAttualeAutoma){
              printf("Innesco (%lu,%lu)\n",statoAttualeAutoma,transizione.second.first);
              
              //Prendo guardia e reset
              auto guardia_reset=this->guardie[transizione.first];
              bool verificaGuardia=guardia_reset.first ? guardia_reset.first(t,statoInizialeODE) : true;
              if(verificaGuardia){
                printf("Guardia!\n");
                
                //Se c'è il reset lo applico
                if(guardia_reset.second) guardia_reset.second(t,statoInizialeODE,statoAttualeAutoma,statoInizialeODE);
                
                statoAttualeAutoma = transizione.second.first;
                locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
                bloccante=false;
                break;
              }else{
              }
            }
          }
          if(bloccante) innescoPronto=true;   //Forzo l'uscita per segnalare l'errore
        }else{
          innesco=gsl_matrix_view_vector(statoInizialeODE,n,1);
          innescoPronto=true;
        }
      }
    }
    if(bloccante){
      printf("Errore all'istante %.10lf: il sistema risulta bloccante nella locazione %lu\n", t, statoAttualeAutoma);
      break;
    }
    
    evoluzioneAutoma.push(pair<double,TipoStato>({t,statoAttualeAutoma}));
    
    infoSimulazione.dinamica=locazioneAttuale->f_ODE;
    infoSimulazione.condizione=locazioneAttuale->condizione;
    infoSimulazione.tCondizione=&istanteCondizione;
    infoSimulazione.indiceCondizione=&indiceCondizione;
    infoSimulazione.t0=t;
    infoSimulazione.T=T-t;
    infoSimulazione.h=infoSimulazione.T/(double)(campioniRimanenti-1); //Risolve il problema dell'indicizzazione, fissando il numero di campioni si evitano problemi numerici
                                                                       //nell'indicizzazione, la variazione rispetto a h nominale è trascurabile
    printf("T: %.10lf\t h: %e\n",infoSimulazione.T,infoSimulazione.h-h);
    
    printf("indiceIstante: %ld\n", indiceIstante);
    
    printf("-------Stato ODE-------\n");
    gsl_vector_fprintf(stdout,statoInizialeODE,"%.10lf");
    
    gsl_matrix* soluzione=mODE(&infoSimulazione,&(innesco.matrix));
    
    //Verifico se la simulazione è stata interrotta, se la condizione non è verificata si può essere interrotta per l'arrivo dell'input
    gsl_vector_view ultimoStato=gsl_matrix_column(soluzione,indiceCondizione);
    printf("---------Ultimo Stato---------\n");
    gsl_vector_fprintf(stdout,&(ultimoStato.vector),"%.10lf");
    bool verificaCondizione= locazioneAttuale->condizione ? locazioneAttuale->condizione(istanteCondizione,&(ultimoStato.vector)) : false;
    printf("uscito a indice: %lu\n", indiceCondizione);
    if(verificaCondizione){
      bloccante=true;
      for(auto transizione : this->transizioni){ //Se il ciclo for termina allora non ci sono transizioni abilitate per uscire dalla locazione, il sistema è bloccante
        if(transizione.first.first == statoAttualeAutoma){
          printf("In (%lu,%lu)\n",statoAttualeAutoma,transizione.second.first);
          
          //Prendo guardia e reset
          auto guardia_reset=this->guardie[transizione.first];
          bool verificaGuardia=guardia_reset.first ? guardia_reset.first(istanteCondizione,&(ultimoStato.vector)) : true;
          if(verificaGuardia){
            printf("Guardia!\n");
            
            //L'innesco non può verificare la condizione di uscita
            double istanteUltimoStatoValido = istanteCondizione-h;
            size_t indiceUltimoStatoValido = indiceCondizione-1;
            gsl_vector_view ultimoStatoValido = gsl_matrix_column(soluzione,indiceUltimoStatoValido);
            
            //Se c'è il reset lo applico
            if(guardia_reset.second) guardia_reset.second(istanteUltimoStatoValido,&(ultimoStatoValido.vector),statoAttualeAutoma,statoInizialeODE);
            else gsl_vector_memcpy(statoInizialeODE,&(ultimoStatoValido.vector));
            
            //Cancellazione stato non valido
            gsl_vector_set_zero(&(ultimoStato.vector));
            statoAttualeAutoma = transizione.second.first;
            bloccante=false;
            break;
          }else{
          }
        }
      }
    }else{
      //Sono uscito perchè è scaduto il periodo, non ci sono stati da eliminare
      indiceCondizione++;
      istanteCondizione=t0+T;
      //gsl_vector_memcpy(statoInizialeODE,&(ultimoStato.vector));
    }
    if(bloccante){
      printf("Errore all'istante %.10lf: il sistema risulta bloccante nella locazione %lu\n", istanteCondizione, statoAttualeAutoma);
      break;
    }
    
    //Aggiornamento cursori e matrici
    gsl_matrix_view parteTotale=gsl_matrix_submatrix(evoluzioneODE,0,indiceIstante,n,soluzione->size2);
    gsl_matrix_add(&(parteTotale.matrix),soluzione);
    indiceIstante += indiceCondizione;
    t=istanteCondizione;
    printf("Istante finale: %.10lf\n",t);
  }
  pair<gsl_matrix*,queue<pair<double,TipoStato>>> res=pair<gsl_matrix*,queue<pair<double,TipoStato>>>({evoluzioneODE,evoluzioneAutoma});
  return res;
}
