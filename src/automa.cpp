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

pair<gsl_matrix*,queue<pair<double,TipoStato>>> Automa::Simulazione(double t0, gsl_vector* x0, TipoStato s0, CodaInput inputs, double T, double h,
                                                                    MetodoODE mODE, MetodoInnesco mInnesco, unsigned passi){
  size_t n=x0->size;
  size_t p=passi-1;
  size_t NumeroCampioni=(size_t)floor(T/h)+1;
  
  struct InfoBaseSimulazione infoSimulazione;
  double istanteCondizione=0.0;
  double t=t0;
  double istanteFinaleSimulazione=t0+((double)(NumeroCampioni-1))*h;
  size_t indiceCondizione=0;
  bool innescoPronto=false;
  bool bloccante=false;
  bool statoIniziale=false;
  
  TipoStato statoAttualeAutoma=s0;
  gsl_vector* statoInizialeODE=gsl_vector_calloc(n);
  gsl_vector_memcpy(statoInizialeODE,x0);
  
  queue<pair<double,TipoStato>> evoluzioneAutoma=queue<pair<double,TipoStato>>();
  gsl_matrix*                   evoluzioneODE=gsl_matrix_calloc(n,NumeroCampioni);
  gsl_matrix_view innesco;
  char argomentiMessaggi[100]={};
  FILE* fileLog=fopen("LogSimulazione.txt","w");
  
  //size_t o=0;
  size_t indiceIstante=0;
  while(indiceIstante < NumeroCampioni-1){
    //  ++o;
    /*------------------Verifica Esistenza Stato Nell'automa------------------*/
    auto locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
    if(locazioneAttuale == this->locazioni.cend()){
      sprintf(argomentiMessaggi,"Non esiste lo stato %lu nell'automa",statoAttualeAutoma);
      LogMessaggio(fileLog,2,t,argomentiMessaggi);
      break;
    }
    
    bloccante=false;
    size_t varZenone=0;
    sprintf(argomentiMessaggi,"Locazione: %lu e Campioni rimanenti: %lu",statoAttualeAutoma,NumeroCampioni-indiceIstante);
    LogMessaggio(fileLog,0,t,argomentiMessaggi);
    
    /*------------------Verifica Stato Iniziale ODE nella locazione------------------*/
    LogMessaggio(fileLog,0,t,"Verifica stato iniziale");
    while(not statoIniziale and varZenone < ZENONE_MAX){
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
        ++varZenone;
      }else{
        sprintf(argomentiMessaggi,"Stato iniziale verificato nella locazione %lu",statoAttualeAutoma);
        LogMessaggio(fileLog,0,t,argomentiMessaggi);
        statoIniziale=true;
      }
      /*++varZenone;
      if(varZenone > ZENONE_MAX){
        LogMessaggio(fileLog,1,t,"Il sistema non riesce a verificare lo stato iniziale, potrebbe essere zenoniano. Simulazione terminata");
        break;
      }*/
    }
    if(bloccante){
      sprintf(argomentiMessaggi,"Il sistema risulta bloccante nella locazione %lu. Simulazione terminata",statoAttualeAutoma);
      LogMessaggio(fileLog,2,t,argomentiMessaggi);
      break;
    }
    if(varZenone >= ZENONE_MAX){
      LogMessaggio(fileLog,1,t,"Il sistema non riesce a verificare lo stato iniziale, potrebbe essere zenoniano. Simulazione terminata");
      break;
    }
    
    /*------------------Preparazione Simulazione------------------*/
    //Preparazione innesco
    infoSimulazione.dinamica=locazioneAttuale->f_ODE;
    infoSimulazione.condizione=locazioneAttuale->condizione;
    infoSimulazione.tCondizione=&istanteCondizione;
    infoSimulazione.indiceCondizione=&indiceCondizione;
    infoSimulazione.t0=t;
    infoSimulazione.h=h;
    
    //Se ci sono degli input quantizzo il loro istante e lo metto come istante finale della simulazione
    double istanteInputQuantizzatoH=inputs.empty() ? -1.0 : floor(inputs.front().first/h)*h;
    double istanteFinale=inputs.empty() ? istanteFinaleSimulazione : istanteInputQuantizzatoH;
    bool innescoInterrottoInput=false;
    
    //Lo stato dell'automa in questo istante è valido
    evoluzioneAutoma.push(pair<double,TipoStato>({t,statoAttualeAutoma}));
    
    //Divisione tra fase di innesco e quella di simulazione
    if(innescoPronto or p==0){ //Se p==0 non devo calcolare altro dell'innesco
      LogMessaggio(fileLog,0,t,"Fase Simulazione");
      infoSimulazione.T=istanteFinale-t;
      if(p==0){
        innescoPronto=true;
        innesco=gsl_matrix_view_vector(statoInizialeODE,n,1);
      }
    }else{
      LogMessaggio(fileLog,0,t,"Fase Innesco");
      size_t campioniRimanenti=(size_t)floor((istanteFinale-t)/h)+1;
      innescoInterrottoInput = campioniRimanenti-1 <= p;
      infoSimulazione.T=(double)min(p,campioniRimanenti-1)*h;
    }
    
    //printf("errore h: %e\n", infoSimulazione.h-h);
    /*------------------Simulazione o Calcolo Innesco------------------*/
    gsl_matrix* soluzione= innescoPronto ? mODE(&infoSimulazione,&(innesco.matrix)) : mInnesco(&infoSimulazione,statoInizialeODE);
    
    /*------------------Verifica Della Causa Del Termine Della Simulazione------------------*/
    //Verifico se la simulazione è stata interrotta, se la condizione non è verificata allora è terminato il periodo
    gsl_vector_view ultimoStato=gsl_matrix_column(soluzione,indiceCondizione);
    bool verificaCondizione= locazioneAttuale->condizione ? locazioneAttuale->condizione(istanteCondizione,&(ultimoStato.vector)) : false;
    if(verificaCondizione){
    
      LogMessaggio(fileLog,0,istanteCondizione,"Condizione di uscita verificata");
      bloccante=true;
      for(auto transizione : this->transizioni){ //Se il ciclo for termina allora non ci sono transizioni abilitate per uscire dalla locazione, il sistema è bloccante
        if(transizione.first.first == statoAttualeAutoma){
        
          //Prendo guardia e reset
          auto guardia_reset=this->guardie[transizione.first];
          bool verificaGuardia=guardia_reset.first ? guardia_reset.first(istanteCondizione,&(ultimoStato.vector)) : true;
          if(verificaGuardia){
          
            sprintf(argomentiMessaggi,"Condizione di guardia verificata per transizione (%lu,%ld)-->%lu",statoAttualeAutoma,transizione.first.second,transizione.second.first);
            LogMessaggio(fileLog,0,istanteCondizione,argomentiMessaggi);
            
            //il nuovo innesco non può verificare la condizione di uscita, e sicuramente non è lo stato iniziale
            size_t indiceUltimoStatoValido = indiceCondizione-1;
            double istanteUltimoStatoValido = t0+((double)(indiceIstante+indiceUltimoStatoValido))*h;
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
      //Se ci sono degli input da elaborare, è arrivato il loro istante
      if((not inputs.empty() and innescoPronto) or (not inputs.empty() and not innescoPronto and innescoInterrottoInput)){
        TipoInput InputAttuale=inputs.front().second;
        sprintf(argomentiMessaggi,"Arrivato input %lu nella locazione %lu",InputAttuale,statoAttualeAutoma);
        LogMessaggio(fileLog,0,istanteCondizione,argomentiMessaggi);
        pair<TipoStato,TipoInput> coppiaChiave=pair<TipoStato,TipoInput>({statoAttualeAutoma,InputAttuale});
        
        //Se esiste la transizione per la locazione attuale
        if(this->transizioni.contains(coppiaChiave)){
        
          bloccante=true;
          //Prendo guardia e reset
          auto guardia_reset=this->guardie[coppiaChiave];
          bool verificaGuardia=guardia_reset.first ? guardia_reset.first(istanteCondizione,&(ultimoStato.vector)) : true;
          if(verificaGuardia){
            auto StatoSuccessivo=this->transizioni[coppiaChiave];
            sprintf(argomentiMessaggi,"Condizione di guardia verificata per transizione (%lu,%ld)-->%lu",statoAttualeAutoma,InputAttuale,StatoSuccessivo.first);
            LogMessaggio(fileLog,0,istanteCondizione,argomentiMessaggi);
            
            //Se c'è il reset lo applico
            if(guardia_reset.second){
              guardia_reset.second(istanteCondizione,&(ultimoStato.vector),statoAttualeAutoma,statoInizialeODE);
              sprintf(argomentiMessaggi,"Reset applicato per transizione (%lu,%ld)-->%lu",statoAttualeAutoma,InputAttuale,StatoSuccessivo.first);
              LogMessaggio(fileLog,0,istanteCondizione,argomentiMessaggi);
            }else{
              gsl_vector_memcpy(statoInizialeODE,&(ultimoStato.vector));
            }
            
            //Cancellazione stato non valido
            statoAttualeAutoma = StatoSuccessivo.first;
            inputs.pop();
            bloccante=false;
          }else{
          }
          
        }else{
          sprintf(argomentiMessaggi,"Non esiste la transizione (%lu,%lu)",statoAttualeAutoma,InputAttuale);
          LogMessaggio(fileLog,2,istanteCondizione,argomentiMessaggi);
          break;
        }
      }
      indiceCondizione++;
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
    
    /*------------------Aggiornamento Soluzioni------------------*/
    //Aggiorno la matrice della soluzione solo se l'innesco è stato interrotto o ho completato una simulazione, dunque se innescoPronto==false
    if(not innescoPronto){
      //Aggiornamento cursori e matrici
      gsl_matrix_view parteTotale=gsl_matrix_submatrix(evoluzioneODE,0,indiceIstante,n,soluzione->size2);
      gsl_matrix_add(&(parteTotale.matrix),soluzione);
      indiceIstante += indiceCondizione;
      t=t0+((double)indiceIstante)*h;
      //printf("Istante finale: %.10lf\n",t);
    }
  }
  
  if(indiceIstante >= NumeroCampioni-1) LogMessaggio(fileLog,0,istanteFinaleSimulazione,"Simulazione terminata con successo");
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
