#include "../include/automa.h"

void Sistema1(double t, gsl_vector* y, gsl_vector* dy);
void Sistema2(double t, gsl_vector* y, gsl_vector* dy);
bool Condizione1(double t, gsl_vector* y);

gsl_matrix* metodo_ODE(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);

int main(int argc, char* argv[]){
  Automa A1({0,1});
  A1.AggiungiRegola(0,0,1,1);
  A1.AggiungiRegola(1,1,0,0);
  A1.ImpostaODE(0,Sistema1);
  A1.ImpostaODE(1,Sistema2);
  A1.ImpostaCondizioneLocazione(0,Condizione1);
  
  double y0[]={2.0};
  gsl_vector_view statoIniziale=gsl_vector_view_array(y0,1);
  queue<pair<double,TipoInput>> codaInput=queue<pair<double,TipoInput>>();
  
  gsl_matrix* simulazione=A1.Simulazione(&(statoIniziale.vector),0,0.0,10.0,1e-2,codaInput,metodo_ODE);
  gsl_matrix_fprintf(stdout,simulazione,"%.12lf");
  gsl_matrix_free(simulazione);
  return 0;
}

bool Condizione1(double t, gsl_vector* y){
  return t > 1.0;
}

bool Condizione2(double t, gsl_vector* y){
  return t > 4.0;
}

void Sistema1(double t, gsl_vector* y, gsl_vector* dy){
  gsl_vector_memcpy(dy,y);
  gsl_vector_scale(dy,-2.0);
}

void Sistema2(double t, gsl_vector* y, gsl_vector* dy){
  gsl_vector_memcpy(dy,y);
  gsl_vector_scale(dy,2.0);
}

gsl_matrix* metodo_ODE(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale){
  return EuleroAvanti(infoSimulazione,statoIniziale);
}
