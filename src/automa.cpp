#include "../include/automa.h"
#include <cmath>
#include <string>
#include <cstring>
#include <stdio.h>
#define ZENONE_MAX 100

void LogMessaggio(FILE* log,int tipo,double istante,const char* messaggio);

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
  bool innescoPronto=false,bloccante=false,statoIniziale=false;
  TipoStato statoAttualeAutoma=s0;
  gsl_vector* statoInizialeODE=gsl_vector_calloc(n);
  gsl_vector_memcpy(statoInizialeODE,x0);
  char argomentiMessaggi[100]={};
  FILE* fileLog=fopen("LogSimulazione.txt","w");
  
  //size_t o=0;
  size_t indiceIstante=0;
  for(;indiceIstante < NumeroCampioni-1;){
    //  ++o;
    bloccante=false;
    auto locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
    size_t campioniRimanenti=NumeroCampioni-indiceIstante,varZenone=0;
    
    sprintf(argomentiMessaggi,"Locazione: %lu e Campioni rimanenti: %lu",statoAttualeAutoma,campioniRimanenti);
    LogMessaggio(fileLog,0,t,argomentiMessaggi);
    
    //Finchè lo stato iniziale non è nella locazione giusta avanzo nell'automa
    LogMessaggio(fileLog,0,t,"Verifica stato iniziale");
    while(not statoIniziale){
      bool verificaCondizione= locazioneAttuale->condizione ? locazioneAttuale->condizione(t,statoInizialeODE) : false;
      if(verificaCondizione){
        
        LogMessaggio(fileLog,0,t,"Condizione di uscita verificata");
        bloccante=true;
        for(auto transizione : this->transizioni){ //Se il ciclo for termina allora non ci sono transizioni abilitate per uscire dalla locazione, il sistema è bloccante
          if(transizione.first.first == statoAttualeAutoma){
            
            //Prendo guardia e reset
            auto guardia_reset=this->guardie[transizione.first];
            bool verificaGuardia=guardia_reset.first ? guardia_reset.first(t,statoInizialeODE) : true;
            if(verificaGuardia){
              sprintf(argomentiMessaggi,"Condizione di guardia verificata per transizione (%lu,%ld)-->%lu",statoAttualeAutoma,transizione.first.second,transizione.second.first);
              LogMessaggio(fileLog,0,t,argomentiMessaggi);
              
              //Se c'è il reset lo applico
              if(guardia_reset.second) {
                guardia_reset.second(t,statoInizialeODE,statoAttualeAutoma,statoInizialeODE);
                sprintf(argomentiMessaggi,"Reset applicato per transizione (%lu,%ld)-->%lu",statoAttualeAutoma,transizione.first.second,transizione.second.first);
                LogMessaggio(fileLog,0,t,argomentiMessaggi);
              }
              
              statoAttualeAutoma = transizione.second.first;
              locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
              bloccante=false;
              break;
            }else{
            }
          }
        }
        if(bloccante) statoIniziale=true;   //Forzo l'uscita per segnalare l'errore
      }else{
        sprintf(argomentiMessaggi,"Stato iniziale verificato nella locazione %lu",statoAttualeAutoma);
        LogMessaggio(fileLog,0,t,argomentiMessaggi);
        statoIniziale=true;
      }
      ++varZenone;
      if(varZenone > ZENONE_MAX){
        LogMessaggio(fileLog,1,t,"Il sistema non riesce a verificare lo stato iniziale, potrebbe essere zenoniano. Simulazione terminata");
        break;
      }
    }
    if(bloccante){
      sprintf(argomentiMessaggi,"Il sistema risulta bloccante nella locazione %lu. Simulazione terminata",statoAttualeAutoma);
      LogMessaggio(fileLog,2,t,argomentiMessaggi);
      break;
    }
    if(varZenone > ZENONE_MAX){
      LogMessaggio(fileLog,1,t,"Il sistema non riesce a verificare lo stato iniziale, potrebbe essere zenoniano. Simulazione terminata");
      break;
    }
    
    //Preparazione innesco
    infoSimulazione.dinamica=locazioneAttuale->f_ODE;
    infoSimulazione.condizione=locazioneAttuale->condizione;
    infoSimulazione.tCondizione=&istanteCondizione;
    infoSimulazione.indiceCondizione=&indiceCondizione;
    infoSimulazione.t0=t;
    
    //Lo stato dell'automa in questo istante è valido
    evoluzioneAutoma.push(pair<double,TipoStato>({t,statoAttualeAutoma}));
    
    //Divisione tra fase di innesco e quella di simulazione
    if(innescoPronto or p==0){ //Se p==0 non devo calcolare altro dell'innesco
      LogMessaggio(fileLog,0,t,"Fase Simulazione");
      infoSimulazione.T=T-t;
      //Risolve il problema dell'indicizzazione, fissando il numero di campioni si evitano problemi numerici nell'indicizzazione, la variazione rispetto a h nominale è trascurabile
      infoSimulazione.h=infoSimulazione.T/(double)(campioniRimanenti-1);
      if(p==0){
        innescoPronto=true;
        innesco=gsl_matrix_view_vector(statoInizialeODE,n,1);
      }
    }else{
      LogMessaggio(fileLog,0,t,"Fase Innesco");
      infoSimulazione.T=((double)p)*h;
      infoSimulazione.h=h;
    }
    //printf("T: %.10lf\t errore h: %e\n",infoSimulazione.T,infoSimulazione.h-h);
    
    //printf("indiceIstante: %ld\n", indiceIstante);
    
    //printf("-------Stato ODE-------\n");
    //gsl_vector_fprintf(stdout,statoInizialeODE,"%.10lf");
    
    //printf("innesco pronto: %d\n",innescoPronto);
    gsl_matrix* soluzione= innescoPronto ? mODE(&infoSimulazione,&(innesco.matrix)) : mInnesco(&infoSimulazione,statoInizialeODE);
    
    //Verifico se la simulazione è stata interrotta, se la condizione non è verificata si può essere interrotta per l'arrivo dell'input
    gsl_vector_view ultimoStato=gsl_matrix_column(soluzione,indiceCondizione);
    //printf("---------Ultimo Stato---------\n");
    //gsl_vector_fprintf(stdout,&(ultimoStato.vector),"%.10lf");
    bool verificaCondizione= locazioneAttuale->condizione ? locazioneAttuale->condizione(istanteCondizione,&(ultimoStato.vector)) : false;
    //printf("uscito a indice: %lu\n", indiceCondizione);
    if(verificaCondizione){
      LogMessaggio(fileLog,0,istanteCondizione,"Condizione di uscita verificata");
      //printf("\tCondizione verificata\n");
      bloccante=true;
      for(auto transizione : this->transizioni){ //Se il ciclo for termina allora non ci sono transizioni abilitate per uscire dalla locazione, il sistema è bloccante
        if(transizione.first.first == statoAttualeAutoma){
          //printf("\tPassaggio (%lu,%lu)\n",statoAttualeAutoma,transizione.second.first);
          
          //Prendo guardia e reset
          auto guardia_reset=this->guardie[transizione.first];
          bool verificaGuardia=guardia_reset.first ? guardia_reset.first(istanteCondizione,&(ultimoStato.vector)) : true;
          if(verificaGuardia){
            sprintf(argomentiMessaggi,"Condizione di guardia verificata per transizione (%lu,%ld)-->%lu",statoAttualeAutoma,transizione.first.second,transizione.second.first);
            LogMessaggio(fileLog,0,istanteCondizione,argomentiMessaggi);
            
            //il nuovo innesco non può verificare la condizione di uscita, e sicuramente non è lo stato iniziale
            double istanteUltimoStatoValido = istanteCondizione-h;
            size_t indiceUltimoStatoValido = indiceCondizione-1;
            gsl_vector_view ultimoStatoValido = gsl_matrix_column(soluzione,indiceUltimoStatoValido);
            
            //Se c'è il reset lo applico
            if(guardia_reset.second){
              guardia_reset.second(istanteUltimoStatoValido,&(ultimoStatoValido.vector),statoAttualeAutoma,statoInizialeODE);
              sprintf(argomentiMessaggi,"Reset applicato per transizione (%lu,%ld)-->%lu",statoAttualeAutoma,transizione.first.second,transizione.second.first);
              LogMessaggio(fileLog,0,istanteCondizione,argomentiMessaggi);
            }else{
              gsl_vector_memcpy(statoInizialeODE,&(ultimoStatoValido.vector));
            }
            
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
    }
    if(bloccante){
      sprintf(argomentiMessaggi,"Il sistema risulta bloccante nella locazione %lu. Simulazione terminata",statoAttualeAutoma);
      LogMessaggio(fileLog,2,istanteCondizione,argomentiMessaggi);
      break;
    }else{
      innescoPronto=not innescoPronto; //Se è pronto l'innesco alla prossima iterazione simulo il sistema
      if(not innescoPronto) statoIniziale=false; //Ho finito la simulazione e devo verificare il nuovo stato iniziale dell'innesco
      if(innescoPronto) innesco=gsl_matrix_submatrix(soluzione,0,0,n,p+1);
    }
    
    //Aggiorno la matrice della soluzione solo se l'innesco è stato interrotto o ho completato una simulazione, dunque se innescoPronto==false
    if(not innescoPronto){
      //Aggiornamento cursori e matrici
      gsl_matrix_view parteTotale=gsl_matrix_submatrix(evoluzioneODE,0,indiceIstante,n,soluzione->size2);
      gsl_matrix_add(&(parteTotale.matrix),soluzione);
      indiceIstante += indiceCondizione;
      t=istanteCondizione;
      //printf("Istante finale: %.10lf\n",t);
    }
  }
  
  if(indiceIstante >= NumeroCampioni-1) LogMessaggio(fileLog,0,t,"Simulazione terminata con successo");
  pair<gsl_matrix*,queue<pair<double,TipoStato>>> res=pair<gsl_matrix*,queue<pair<double,TipoStato>>>({evoluzioneODE,evoluzioneAutoma});
  fclose(fileLog);
  return res;
}

void LogMessaggio(FILE* log,int tipo,double istante,const char* messaggio){
  char strTipo[15]={};
  switch(tipo){
    case 0:
      strcpy(strTipo,"[INFO] ");
      break;
    case 1:
      strcpy(strTipo,"[WARNING] ");
      break;
    case 2:
      strcpy(strTipo,"[ERROR] ");
      break;
  }
  
  char strIstante[30]={};
  sprintf(strIstante,"in istante %e: ",istante);
  string messaggioLog=string();
  messaggioLog += strTipo;
  messaggioLog += strIstante;
  messaggioLog += messaggio;
  
  fprintf(log,"%s\n",messaggioLog.c_str());
}
