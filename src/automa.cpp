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
            gsl_matrix* (*metodo_ODE)(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale)){

  const size_t n=y0->size;
  TipoStato statoAttualeAutoma=s0;
  auto locazioneAttuale=this->locazioni.find(s0);
  gsl_vector* statoInizialeSimulazione=gsl_vector_calloc(n);
  double t=t0,istanteCondizione=0.0;
  size_t indiceIstante=0;
  size_t NumeroCampioni=(size_t)floor(T/h)+1,indiceCondizione=0;
  gsl_vector_memcpy(statoInizialeSimulazione,y0);
  
  //Matrice del calcolo totale
  gsl_matrix* O_sim_totale=gsl_matrix_calloc(n,NumeroCampioni);
  
  //Simulazione a blocchi del sistema ibrido
  struct InfoBaseSimulazione infoSimulazione;
  for(;indiceIstante+1 < NumeroCampioni;indiceIstante += indiceCondizione+1){
    //Impostazione della simulazione con la locazione attuale
    infoSimulazione.dinamica=locazioneAttuale->f_ODE;
    infoSimulazione.condizione=locazioneAttuale->condizione;
    infoSimulazione.tCondizione=&istanteCondizione;
    infoSimulazione.indiceCondizione=&indiceCondizione;
    infoSimulazione.t0=t;
    infoSimulazione.T= seqInput.empty() ? T-(t-t0) : seqInput.front().first-(t-t0);
    infoSimulazione.h=h;
    
    gsl_matrix* O_sim=metodo_ODE(&infoSimulazione,statoInizialeSimulazione);
    if(O_sim->size2 > 0){
      gsl_vector_view ultimoStato=gsl_matrix_column(O_sim,indiceCondizione);
      gsl_matrix_view parte_O_sim_totale=gsl_matrix_submatrix(O_sim_totale,0,indiceIstante,n,O_sim->size2);
      gsl_matrix_add(&(parte_O_sim_totale.matrix),O_sim);
      
      //Vedo se sono uscito per la condiizone, mi serve l'ultimo stato
      bool verificaCondizione=locazioneAttuale->condizione==nullptr ? false : locazioneAttuale->condizione(istanteCondizione,&(ultimoStato.vector));
      if(verificaCondizione){
        //È verificata la condizione di uscita, prendo la prima transizione disponibile
        for(auto transizione : this->transizioni){
          if(transizione.first.first == statoAttualeAutoma){
            //prendo la sua condizione di guardia e la verifico
            auto guardia_reset=this->guardie[transizione.first];
            bool verificaGuardia = guardia_reset.first==nullptr ? true : guardia_reset.first(istanteCondizione,&(ultimoStato.vector));
            if(verificaGuardia){
              //Transizione dell'automa
              statoAttualeAutoma=transizione.second.first;
              locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
              
              //Se c'è il reset lo eseguo
              if(guardia_reset.second) guardia_reset.second(istanteCondizione,&(ultimoStato.vector),statoInizialeSimulazione);
              else gsl_vector_memcpy(statoInizialeSimulazione,&(ultimoStato.vector));
              break;
            }
          }
        }
      }else if(not seqInput.empty()){
        auto Istante_Input=seqInput.front();
        pair<TipoStato,TipoInput> coppia=pair<TipoStato,TipoInput>({statoAttualeAutoma,Istante_Input.second});
        pair<TipoStato,TipoOutput> transizione=this->transizioni[coppia];
        //prendo la sua condizione di guardia e la verifico
        auto guardia_reset=this->guardie[coppia];
        bool verificaGuardia = guardia_reset.first==nullptr ? true : guardia_reset.first(Istante_Input.first,&(ultimoStato.vector));
        if(verificaGuardia){
          //Transizione dell'automa
          statoAttualeAutoma=transizione.first;
          locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
          
          //Se c'è il reset lo eseguo
          if(guardia_reset.second) guardia_reset.second(Istante_Input.first,&(ultimoStato.vector),statoInizialeSimulazione);
          else gsl_vector_memcpy(statoInizialeSimulazione,&(ultimoStato.vector));
          seqInput.pop();
        }
      }
      t = istanteCondizione+h;
      gsl_matrix_free(O_sim);
    }
  }
  return O_sim_totale;
}
