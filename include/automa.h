#include <set>
#include <map>
#include <queue>
#include <gsl/gsl_matrix.h>

using namespace std;

typedef unsigned long TipoStato;
typedef long TipoInput;
typedef long TipoOutput;

class Locazione{
  TipoStato stato;
  void (*f_ODE)(double,gsl_vector*,gsl_vector*);
  bool (*condizione)(double,gsl_vector*);
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
  map<pair<TipoStato,TipoInput>,pair<bool (*)(double,gsl_vector*),void (*)(double,gsl_vector*,gsl_vector*)>> guardie;
public:
  Automa(set<TipoStato> stati);
  Automa(set<TipoStato>& stati);
  bool AggiungiRegola(TipoStato stato, TipoInput input, TipoStato statoSucc, TipoOutput output);
  pair<queue<TipoStato>,queue<TipoOutput>> ValutaInput(queue<TipoInput> stringaInput,TipoStato s0);
  void ImpostaODE(TipoStato s, void (*f)(double,gsl_vector*,gsl_vector*));
  void ImpostaCondizioneLocazione(TipoStato s, bool (*f)(double,gsl_vector*));
  void ImpostaCondizioneGuardia(TipoStato s, TipoInput i, bool (*g)(double,gsl_vector*),void (*r)(double,gsl_vector*,gsl_vector*));
  
  pair<gsl_matrix*,queue<pair<double,TipoStato>>> Simulazione(gsl_vector* y0, TipoStato s0, double t0, double T, double h, queue<pair<double,TipoInput>> seqInput,
  gsl_matrix* (*metodo_ODE)(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, gsl_vector* y0,bool (*condizione)(double,gsl_vector*), double *tCondizione,unsigned* iCondizione));
};
