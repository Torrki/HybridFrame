#include "../include/automa.h"
#include <ode.h>
#include <cmath>

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
    this->guardie[chiave]=pair<bool (*)(double,gsl_vector*),void (*)(double,gsl_vector*,gsl_vector*)>({nullptr,nullptr});
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

void Automa::ImpostaODE(TipoStato s, void (*f)(double,gsl_vector*,gsl_vector*)){
  if(this->locazioni.contains(s)){
    auto NodoLocazione=this->locazioni.extract(s);
    NodoLocazione.value().f_ODE=f;
    this->locazioni.insert(move(NodoLocazione));
  }
}

void Automa::ImpostaCondizioneLocazione(TipoStato s, bool (*f)(double,gsl_vector*)){
  if(this->locazioni.contains(s)){
    auto NodoLocazione=this->locazioni.extract(s);
    NodoLocazione.value().condizione=f;
    this->locazioni.insert(move(NodoLocazione));
  }
}

void Automa::ImpostaCondizioneGuardia(TipoStato s, TipoInput i, bool (*g)(double,gsl_vector*),void (*r)(double,gsl_vector*,gsl_vector*)){
  pair<TipoStato,TipoInput> coppia({s,i});
  if(this->transizioni.contains(coppia)){
    this->guardie[coppia]=pair<bool (*)(double,gsl_vector*),void (*)(double,gsl_vector*,gsl_vector*)>({g,r});
  }
}

gsl_matrix* Automa::Simulazione(gsl_vector* y0, TipoStato s0, double t0, double T, double h, queue<pair<double,TipoInput>> seqInput,
            gsl_matrix* (*metodo_ODE)(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, gsl_vector* y0,bool (*condizione)(double,gsl_vector*), double *tCondizione)){

  TipoStato statoAttuale=s0;
  auto locazioneAttuale=this->locazioni.find(s0);
  double t=t0;
  unsigned indiceIstante=0;
  unsigned NCampioni=(unsigned)floor(T/h)+1;
  
  gsl_vector* statoInizialeSim=gsl_vector_alloc(y0->size);
  
  
  gsl_vector_memcpy(statoInizialeSim,y0);
  gsl_matrix* O_sim_totale=gsl_matrix_calloc(y0->size,NCampioni);
  
  while(t < t0+T){
    //Ottengo i parametri per il solver
    auto f=locazioneAttuale->f_ODE;
    auto condizioneLoc=locazioneAttuale->condizione;
    double ultimoIstante=0.0;
    unsigned ultimoIndice=0;
    
    //Prendo l'istante dell'input come istante finale per la simulazione
    pair<double,TipoInput> istanteInput = seqInput.empty() ? pair<double,TipoInput>({-1,0}) : seqInput.front();
    double istanteFinaleSimulazione = seqInput.empty() ? t0+T : floor(istanteInput.first/h)*h;  //Prendo l'istante più vicino all'istante dell'input per rientrare nella griglia di campionamento
    
    gsl_matrix* O_sim=metodo_ODE(f,t,istanteFinaleSimulazione-t,statoInizialeSim,condizioneLoc,&ultimoIstante);
    
    //Aggiungo alla soluzione totale
    ultimoIndice=floor((ultimoIstante-t)/h);
    gsl_matrix_view subMat_totale=gsl_matrix_submatrix(O_sim_totale,0,indiceIstante,y0->size,min(ultimoIndice+1,NCampioni-indiceIstante));
    gsl_matrix_view subMat_sim=gsl_matrix_submatrix(O_sim,0,0,y0->size,min(ultimoIndice+1,NCampioni-indiceIstante));
    gsl_matrix_add(&(subMat_totale.matrix),&(subMat_sim.matrix));
    
    t = ultimoIstante;
    indiceIstante += ultimoIndice+1;
    gsl_vector_view ultimoStato=gsl_matrix_column(O_sim,ultimoIndice);
    
    //Verifico se sono uscito per la condizione
    bool verificaCondizione = condizioneLoc == nullptr ? false : condizioneLoc(t,&(ultimoStato.vector));
    if(verificaCondizione){
      //La simulazione deve riprendere da ultimoIstante+1
      
      //Prendo la prima transizione abilitata nella mappa
      for(auto transizione : this->transizioni){
        if(transizione.first.first == statoAttuale){
          //Verifico la condizione di guardia se esiste
          auto guardia_reset=this->guardie[transizione.first];
          bool transizioneAbilitata= guardia_reset.first == nullptr ? true : guardia_reset.first(t,&(ultimoStato.vector));
          if(transizioneAbilitata){
            //Nuovo stato
            statoAttuale=transizione.second.first;
            locazioneAttuale=this->locazioni.find(statoAttuale);
            
            //Prendo il reset
            if(guardia_reset.second == nullptr) gsl_vector_memcpy(statoInizialeSim,&(ultimoStato.vector));
            else guardia_reset.second(t,&(ultimoStato.vector),statoInizialeSim);
            break;
          }
        }
      }
    }else if(not seqInput.empty()){
      //Verifico la condizione di guardia se esiste
      pair<TipoStato,TipoInput> chiave({statoAttuale,istanteInput.second});
      auto guardia_reset=this->guardie[chiave];
      bool transizioneAbilitata= guardia_reset.first == nullptr ? true : guardia_reset.first(t,&(ultimoStato.vector));
      if(transizioneAbilitata){
        //Nuovo stato
        statoAttuale=this->transizioni[chiave].first;
        locazioneAttuale=this->locazioni.find(statoAttuale);
        
        //Prendo il reset
        if(guardia_reset.second == nullptr) gsl_vector_memcpy(statoInizialeSim,&(ultimoStato.vector));
        else guardia_reset.second(t,&(ultimoStato.vector),statoInizialeSim);
        seqInput.pop();
      }
    }
    gsl_matrix_free(O_sim);
  }
  gsl_vector_free(statoInizialeSim);
  
  return O_sim_totale;
}
