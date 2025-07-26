#include <set>
#include <map>
#include <queue>
#include <ode.h>
#include <gsl/gsl_matrix.h>

using namespace std;

typedef unsigned long TipoStato;
typedef long TipoInput;
typedef long TipoOutput;

class Locazione{
  TipoStato stato;
  ODE f_ODE;
  Condizione condizione;
public:
  Locazione(TipoStato s);
  TipoStato GetStato(){return this->stato;}
  
  //Operatore per usare il tipo set con Locazione
  bool operator<(const Locazione& altro) const{return this->stato < altro.stato;}
  friend class Automa;
};

class Automa{
  set<Locazione> locazioni;
  map<pair<TipoStato,TipoInput>,pair<TipoStato,TipoOutput>> transizioni;
  map<pair<TipoStato,TipoInput>,pair<Condizione,void (*)(double,gsl_vector*,gsl_vector*)>> guardie;
public:
  Automa(set<TipoStato> stati);
  Automa(set<TipoStato>& stati);
  bool AggiungiRegola(TipoStato stato, TipoInput input, TipoStato statoSucc, TipoOutput output);
  pair<queue<TipoStato>,queue<TipoOutput>> ValutaInput(queue<TipoInput> stringaInput,TipoStato s0);
  void ImpostaODE(TipoStato s, ODE f);
  void ImpostaCondizioneLocazione(TipoStato s, Condizione f);
  void ImpostaCondizioneGuardia(TipoStato s, TipoInput i, Condizione g,void (*r)(double,gsl_vector*,gsl_vector*));
  
  gsl_matrix* Simulazione(gsl_vector* y0, TipoStato s0, double t0, double T, double h, queue<pair<double,TipoInput>> seqInput,
              gsl_matrix* (*metodo_ODE)(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale));
};
