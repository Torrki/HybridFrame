#include <set>
#include <map>
#include <queue>
#include <ode.h>
#include <gsl/gsl_matrix.h>

using namespace std;

typedef unsigned long TipoStato;
typedef long TipoInput;
typedef long TipoOutput;
typedef gsl_matrix* (*MetodoODE)(struct InfoBaseSimulazione*,gsl_matrix*);
typedef gsl_matrix* (*MetodoInnesco)(struct InfoBaseSimulazione*,gsl_vector*);
typedef void (*FunzioneReset)(double,gsl_vector*,TipoStato,gsl_vector*);

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
  map<pair<TipoStato,TipoInput>,pair<Condizione,FunzioneReset>> guardie;
public:
  Automa(set<TipoStato> stati);
  Automa(set<TipoStato>& stati);
  bool AggiungiRegola(TipoStato stato, TipoInput input, TipoStato statoSucc, TipoOutput output);
  pair<queue<TipoStato>,queue<TipoOutput>> ValutaInput(queue<TipoInput> stringaInput,TipoStato s0);
  void ImpostaODE(TipoStato s, ODE f);
  void ImpostaCondizioneLocazione(TipoStato s, Condizione f);
  void ImpostaCondizioneGuardia(TipoStato s, TipoInput i, Condizione g,FunzioneReset r);
  pair<gsl_matrix*,queue<pair<double,TipoStato>>> Simulazione(double t0, gsl_vector* x0, TipoStato s0, double T, double h, MetodoODE ODE, MetodoInnesco Innesco, unsigned passi);
};
