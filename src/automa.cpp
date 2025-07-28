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

pair<gsl_matrix*,queue<pair<double,TipoStato>>> Automa::Simulazione(gsl_vector* y0, TipoStato s0, double t0, double T, double h, queue<pair<double,TipoInput>> seqInput,
            gsl_matrix* (*metodo_ODE)(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale)){

  const double ERRORE_PASSO_RIFERIMENTO=erf(5e-4);
  if(erf(h) < ERRORE_PASSO_RIFERIMENTO){
    printf("Warning: h è troppo piccolo, molto probabili errori macchina dovuta alla poca precisione\n");
  }
  const size_t n=y0->size;
  TipoStato statoAttualeAutoma=s0;
  auto locazioneAttuale=this->locazioni.find(s0);
  gsl_vector* statoInizialeSimulazione=gsl_vector_calloc(n);
  double t=t0,istanteCondizione=0.0;
  size_t indiceIstante=0;
  double divisioneLog=log10(T)-log10(h);
  size_t NumeroCampioni=(size_t)ceil(pow(10.0,divisioneLog)),indiceCondizione=0;
  //printf("Numero Campioni: %lu\n",NumeroCampioni);
  gsl_vector_memcpy(statoInizialeSimulazione,y0);
  //printf("erf: %.12lf\n",erf(h));
  
  //Matrice del calcolo totale
  gsl_matrix* O_sim_totale=gsl_matrix_calloc(n,NumeroCampioni);
  queue<pair<double,TipoStato>> evoluzioneAutoma=queue<pair<double,TipoStato>>();
  evoluzioneAutoma.push(pair<double,TipoStato>({t0,s0}));
  
  //Simulazione a blocchi del sistema ibrido
  struct InfoBaseSimulazione infoSimulazione;
  for(;indiceIstante < NumeroCampioni and T-t-t0 > h;indiceIstante += indiceCondizione){
  unsigned p=0;
  //for(;p < 20 and T-t-t0 > h;indiceIstante += indiceCondizione){
    ++p;
    double istanteFinaleSimulazione=seqInput.empty() ? T+t0 : round(seqInput.front().first/h)*h; //Se ci sono input prendo l'istante di tempo nella griglia più vicino a quello nominale
    //Impostazione della simulazione con la locazione attuale
    infoSimulazione.dinamica=locazioneAttuale->f_ODE;
    infoSimulazione.condizione=locazioneAttuale->condizione;
    infoSimulazione.tCondizione=&istanteCondizione;
    infoSimulazione.indiceCondizione=&indiceCondizione;
    infoSimulazione.t0=t;
    infoSimulazione.T= istanteFinaleSimulazione-t;
    infoSimulazione.h=h;
    
    gsl_matrix* O_sim=metodo_ODE(&infoSimulazione,statoInizialeSimulazione);
    gsl_matrix_view parte_O_sim_totale=gsl_matrix_submatrix(O_sim_totale,0,indiceIstante,n,O_sim->size2);
    //printf("Sub\n");
    gsl_matrix_add(&(parte_O_sim_totale.matrix),O_sim);
    gsl_vector_view ultimoStato=gsl_matrix_column(&(parte_O_sim_totale.matrix),indiceCondizione);
    
    //Vedo se sono uscito per la condiizone, mi serve l'ultimo stato
    bool verificaCondizione=locazioneAttuale->condizione==nullptr ? false : locazioneAttuale->condizione(istanteCondizione,&(ultimoStato.vector));
    //printf("Simulato\tindiceCondizione:%lu\tstatoAutoma:%lu\tcondizione:%d\tistanteCondizione: %.10lf\tdiff: %d\n",indiceCondizione,statoAttualeAutoma,verificaCondizione,istanteCondizione,T-t-t0 < h);
    if(verificaCondizione){
      //È verificata la condizione di uscita, prendo la prima transizione disponibile
      for(auto transizione : this->transizioni){
        if(transizione.first.first == statoAttualeAutoma){
          //prendo la sua condizione di guardia e la verifico
          auto guardia_reset=this->guardie[transizione.first];
          bool verificaGuardia = guardia_reset.first==nullptr ? true : guardia_reset.first(istanteCondizione,&(ultimoStato.vector));
          //printf("Guardia: %d\n",verificaGuardia);
          if(verificaGuardia){
            //Transizione dell'automa
            statoAttualeAutoma=transizione.second.first;
            locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
            
            //l'ultimo stato non soddisfa la condizione, non può stare nella soluzione numerica
            //Se c'è il reset lo eseguo
            unsigned indiceUltimoStatoValido= indiceCondizione == 0 ? 0 : indiceCondizione-1;
            //se indiceCondizione è 0 allora già alla partenza lo stato dell'automa non è valido
            if(indiceCondizione==0) evoluzioneAutoma.pop();
            double istanteUltimoStatoValido=fmax(infoSimulazione.t0,istanteCondizione-h);
            gsl_vector_view ultimoStatoValido=gsl_matrix_column(&(parte_O_sim_totale.matrix),indiceUltimoStatoValido);
            //gsl_vector_fprintf(stdout,&(ultimoStato.vector),"%.10lf");
            if(guardia_reset.second) guardia_reset.second(istanteUltimoStatoValido,&(ultimoStatoValido.vector),statoInizialeSimulazione);
            else gsl_vector_memcpy(statoInizialeSimulazione,&(ultimoStatoValido.vector));
            gsl_vector_set_all(&(ultimoStato.vector),0.0);
            gsl_vector_set_all(&(ultimoStatoValido.vector),0.0);
            istanteCondizione=istanteUltimoStatoValido;
            indiceCondizione=indiceUltimoStatoValido;
            evoluzioneAutoma.push(pair<double,TipoStato>({istanteCondizione,statoAttualeAutoma}));
            break;
          }
        }
      }
    }else if(not seqInput.empty()){
      auto Istante_Input=seqInput.front();
      pair<TipoStato,TipoInput> coppia=pair<TipoStato,TipoInput>({statoAttualeAutoma,Istante_Input.second});
      if(this->transizioni.contains(coppia)){ //Se esiste la transizione statoAttuale-input
        //printf("Input\n");
        pair<TipoStato,TipoOutput> transizione=this->transizioni[coppia];
        //prendo la sua condizione di guardia e la verifico
        auto guardia_reset=this->guardie[coppia];
        bool verificaGuardia = guardia_reset.first==nullptr ? true : guardia_reset.first(Istante_Input.first,&(ultimoStato.vector));
        if(verificaGuardia){
          //Transizione dell'automa
          statoAttualeAutoma=transizione.first;
          locazioneAttuale=this->locazioni.find(statoAttualeAutoma);
          
          //l'ultimo stato è quello dell'istante dell'input, può stare nella soluzione numerica
          //unsigned indiceUltimoStatoValido= indiceCondizione == 0 ? 0 : indiceCondizione-1;
          //double istanteUltimoStatoValido=fmax(infoSimulazione.t0,istanteCondizione-h);
          //gsl_vector_view ultimoStatoValido=gsl_matrix_column(O_sim,indiceUltimoStatoValido);
          if(guardia_reset.second) guardia_reset.second(istanteCondizione,&(ultimoStato.vector),statoInizialeSimulazione);
          else gsl_vector_memcpy(statoInizialeSimulazione,&(ultimoStato.vector));
          gsl_vector_set_all(&(ultimoStato.vector),0.0);
          seqInput.pop();
          evoluzioneAutoma.push(pair<double,TipoStato>({istanteCondizione,statoAttualeAutoma}));
          //++indiceCondizione;
        }
      }else{
        fprintf(stderr,"Errore all'istante %.12lf: Non esiste la transizione (%lu,%lu)\n", Istante_Input.first, statoAttualeAutoma,Istante_Input.second);
        break;
      }
    }
    t = istanteCondizione;
    gsl_matrix_free(O_sim);
  }
  return pair<gsl_matrix*,queue<pair<double,TipoStato>>>({O_sim_totale,evoluzioneAutoma});
}
