#include "../include/automa.h"

void Sistema1(double t, gsl_vector* y, gsl_vector* dy);
void Sistema2(double t, gsl_vector* y, gsl_vector* dy);
bool Condizione1(double t, gsl_vector* y);
bool Condizione2(double t, gsl_vector* y);

gsl_matrix* metodo_ODE(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);

struct InfoSimulazione{
  size_t dimensioneInfo;
  size_t righe,colonne,numeroStatiDiscreti;
  double t0,T;
};

int main(int argc, char* argv[]){
  Automa A1({0,1});
  A1.AggiungiRegola(0,0,1,1);
  A1.AggiungiRegola(1,1,0,0);
  A1.ImpostaODE(0,Sistema1);
  A1.ImpostaODE(1,Sistema2);
  A1.ImpostaCondizioneLocazione(0,Condizione1);
  A1.ImpostaCondizioneLocazione(1,Condizione2);
  
  double y0[]={5.0};
  gsl_vector_view statoIniziale=gsl_vector_view_array(y0,1);
  queue<pair<double,TipoInput>> codaInput=queue<pair<double,TipoInput>>();
  codaInput.push(pair<double,TipoInput>({1.0,1}));
  codaInput.push(pair<double,TipoInput>({1.5,1}));
  
  TipoStato s0=1;
  double t0=0.0,T=2.0,h=1e-6;
  auto risultatoSimulazione=A1.Simulazione(&(statoIniziale.vector),s0,t0,T,h,codaInput,metodo_ODE);
  gsl_matrix* simulazione=risultatoSimulazione.first;
  gsl_matrix* evoluzioneAutoma=gsl_matrix_alloc(2,risultatoSimulazione.second.size());
  struct InfoSimulazione info={.dimensioneInfo=sizeof(struct InfoSimulazione),.righe=simulazione->size1,.colonne=simulazione->size2,.numeroStatiDiscreti=risultatoSimulazione.second.size(),.t0=t0,.T=T};
  
  for(unsigned k=0;not risultatoSimulazione.second.empty();++k){
    auto coppia=risultatoSimulazione.second.front();
    gsl_matrix_set(evoluzioneAutoma,0,k,coppia.first);
    gsl_matrix_set(evoluzioneAutoma,1,k,coppia.second);
    //printf("istante: %.12lf\tstato: %lu\n",coppia.first,coppia.second);
    risultatoSimulazione.second.pop();
  }
  
  FILE* fileSimulazione=fopen("datiSimulazione","wb");
  fwrite(&info,sizeof(struct InfoSimulazione),1,fileSimulazione);
  gsl_matrix_fwrite(fileSimulazione,simulazione);
  gsl_matrix_fwrite(fileSimulazione,evoluzioneAutoma);
  fclose(fileSimulazione);
  
  gsl_matrix_free(simulazione);
  gsl_matrix_free(evoluzioneAutoma);
  return 0;
}

bool Condizione1(double t, gsl_vector* y){
  double vc=gsl_vector_get(y,0);
  return vc < 1.0;
}

bool Condizione2(double t, gsl_vector* y){
  double vc=gsl_vector_get(y,0);
  return vc > 4.0;
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
  //printf("Numero Campioni: %u\n",(unsigned)floor(infoSimulazione->T/infoSimulazione->h)+1);
  unsigned stadi=4;
  double A[stadi*stadi]={0.0,0.0,0.0,0.0, 5e-1,0.0,0.0,0.0, 0.0,5e-1,0.0,0.0, 0.0,0.0,1.0,0.0};
  double B[stadi]={1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
  return RungeKuttaEsplicito(infoSimulazione,A,B,stadi,statoIniziale);
}
