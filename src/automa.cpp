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

pair<gsl_matrix*,queue<pair<double,TipoStato>>> Automa::Simulazione(gsl_vector* y0, TipoStato s0, double t0, double T, double h, queue<pair<double,TipoInput>> seqInput,
            gsl_matrix* (*metodo_ODE)(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, gsl_vector* y0,bool (*condizione)(double,gsl_vector*), double *tCondizione, unsigned* iCondizione)){

  //Impostazione simulazione
  TipoStato statoAttuale=s0;
  auto locazioneAttuale=this->locazioni.find(s0);
  double t=t0;
  unsigned indiceSimulazione=0;
  unsigned NumeroCampioni=(unsigned)floor(T/h)+1; //T è il periodo di simulazione, t0+T è l'istante finale
  gsl_vector* statoInizialeSim=gsl_vector_alloc(y0->size);
  gsl_vector_memcpy(statoInizialeSim,y0);
  gsl_matrix* O_sim_totale=gsl_matrix_calloc(y0->size,NumeroCampioni);
  queue<pair<double,TipoStato>> codaStatiDiscreti=queue<pair<double,TipoStato>>();
  codaStatiDiscreti.push(pair<double,TipoStato>({t0,s0}));
  
  while(t < t0+T){
    //Ottengo i parametri per il solver
    auto f=locazioneAttuale->f_ODE;
    auto condizioneLoc=locazioneAttuale->condizione;
    double ultimoIstante=0.0;
    unsigned ultimoIndice=0;
    
    //Prendo l'istante dell'input come istante finale per la simulazione
    pair<double,TipoInput> Istante_Input = seqInput.empty() ? pair<double,TipoInput>({-1.0,0}) : seqInput.front();
    //Prendo l'istante più vicino all'istante dell'input per rientrare nella griglia di campionamento
    double istanteFinaleSimulazione = seqInput.empty() ? t0+T : floor(Istante_Input.first/h)*h;
    
    gsl_matrix* O_sim=metodo_ODE(f,t,istanteFinaleSimulazione-t,statoInizialeSim,condizioneLoc,&ultimoIstante,&ultimoIndice);
    
    //Aggiungo alla soluzione totale
    //ultimoIndice=floor((ultimoIstante-t)/h);
    gsl_matrix_view subMat_totale=gsl_matrix_submatrix(O_sim_totale,0,indiceSimulazione,y0->size,ultimoIndice+1);
    gsl_matrix_view subMat_sim=gsl_matrix_submatrix(O_sim,0,0,y0->size,ultimoIndice+1);
    gsl_matrix_add(&(subMat_totale.matrix),&(subMat_sim.matrix));
    t=ultimoIstante;
    indiceSimulazione += ultimoIndice;
    
    gsl_vector_view ultimoStato=gsl_matrix_column(O_sim_totale,indiceSimulazione);
    
    //Verifico se sono uscito per la condizione
    bool verificaCondizione = condizioneLoc == nullptr ? false : condizioneLoc(ultimoIstante,&(ultimoStato.vector));
    if(verificaCondizione){
      //La simulazione deve riprendere da ultimoIndice+1
      
      //Prendo la prima transizione abilitata nella mappa
      for(auto transizione : this->transizioni){
        if(transizione.first.first == statoAttuale){
          //Verifico la condizione di guardia se esiste
          auto guardia_reset=this->guardie[transizione.first];
          bool transizioneAbilitata= guardia_reset.first == nullptr ? true : guardia_reset.first(ultimoIstante,&(ultimoStato.vector));
          if(transizioneAbilitata){
            //Nuovo stato
            statoAttuale=transizione.second.first;
            locazioneAttuale=this->locazioni.find(statoAttuale);
            codaStatiDiscreti.push(pair<double,TipoStato>({t,statoAttuale}));
            
            //Prendo il reset, ultimoIstante e ultimoStato verificano la condizione di uscita dalla locazione
            //ultimoStato non può stare nella soluzione numerica, uso il punultimo stato, l'ultimo che non verifica la condizione, quello in t-2h
            gsl_vector_set_all(&(ultimoStato.vector),0.0);
            ultimoStato= indiceSimulazione == 0 ? ultimoStato : gsl_matrix_column(O_sim_totale,indiceSimulazione-1);
            if(guardia_reset.second == nullptr) gsl_vector_memcpy(statoInizialeSim,&(ultimoStato.vector));
            else guardia_reset.second(ultimoIstante-h,&(ultimoStato.vector),statoInizialeSim);
            
            //Aggiornamento per la nuova simulazione, dall'istante successivo a quello di uscita
            //t += h;
            //indiceSimulazione += 1;
            break;
          }
        }
      }
    }else if(not seqInput.empty()){
      //Verifico la condizione di guardia se esiste
      pair<TipoStato,TipoInput> chiave({statoAttuale,Istante_Input.second});
      auto guardia_reset=this->guardie[chiave];
      bool transizioneAbilitata= guardia_reset.first == nullptr ? true : guardia_reset.first(ultimoIstante,&(ultimoStato.vector));
      if(transizioneAbilitata){
        //Nuovo stato
        statoAttuale=this->transizioni[chiave].first;
        locazioneAttuale=this->locazioni.find(statoAttuale);
        codaStatiDiscreti.push(pair<double,TipoStato>({t,statoAttuale}));
        
        //Prendo il reset
        gsl_vector_set_all(&(ultimoStato.vector),0.0);
        ultimoStato= indiceSimulazione == 0 ? ultimoStato : gsl_matrix_column(O_sim_totale,indiceSimulazione-1);
        if(guardia_reset.second == nullptr) gsl_vector_memcpy(statoInizialeSim,&(ultimoStato.vector));
        else guardia_reset.second(ultimoIstante-h,&(ultimoStato.vector),statoInizialeSim);
        seqInput.pop();
        
        //Aggiornamento per la nuova simulazione, dall'istante successivo a quello di uscita
        //t += h;
        //indiceSimulazione += 1;
      }
    }
    gsl_matrix_free(O_sim);
  }
  gsl_vector_free(statoInizialeSim);
  
  return pair<gsl_matrix*,queue<pair<double,TipoStato>>>({O_sim_totale,codaStatiDiscreti});
}
