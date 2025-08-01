#include "../include/automa.h"

void Sistema1(double t,gsl_vector* y,gsl_vector* dy);
void Sistema2(double t,gsl_vector* y,gsl_vector* dy);
void Sistema3(double t,gsl_vector* y,gsl_vector* dy);
bool Condizione1(double t, gsl_vector* y);
bool Condizione2(double t, gsl_vector* y);
bool Condizione3(double t, gsl_vector* y);

gsl_matrix* metodoPerODE(struct InfoBaseSimulazione* info,gsl_matrix* innesco);
gsl_matrix* metodoPerInnesco(struct InfoBaseSimulazione* info,gsl_vector* x0);

struct InfoSimulazione{
  size_t dimensioneInfo;
  size_t righe,colonne,numeroStatiDiscreti;
  double t0,T;
};

int main(int argc, char* argv[]){
  Automa SistemaIbrido=Automa({0,1,2});
  SistemaIbrido.ImpostaODE(0,Sistema1);
  SistemaIbrido.ImpostaODE(1,Sistema2);
  SistemaIbrido.ImpostaODE(2,Sistema3);
  SistemaIbrido.ImpostaCondizioneLocazione(0,Condizione1);
  SistemaIbrido.ImpostaCondizioneLocazione(1,Condizione2);
  SistemaIbrido.ImpostaCondizioneLocazione(2,Condizione3);
  
  SistemaIbrido.AggiungiRegola(0,2,2,0);
  SistemaIbrido.AggiungiRegola(2,1,1,2);
  SistemaIbrido.AggiungiRegola(1,2,2,1);
  
  double x0[]={5e-1,-5e-1};
  TipoStato s0=0;
  double t0=0.0,T=2.1,h=1e-3;
  
  gsl_vector_view x0_vect=gsl_vector_view_array(x0,2);
  auto risultato=SistemaIbrido.Simulazione(t0,&(x0_vect.vector),s0,T,h,metodoPerODE,metodoPerInnesco,1);
  
  gsl_matrix* simulazione=risultato.first;
  gsl_matrix* evoluzioneAutoma=gsl_matrix_alloc(2,risultato.second.size());
  struct InfoSimulazione info={.dimensioneInfo=sizeof(struct InfoSimulazione),.righe=simulazione->size1,
                              .colonne=simulazione->size2,.numeroStatiDiscreti=risultato.second.size(),.t0=t0,.T=T};
  
  for(unsigned k=0;not risultato.second.empty();++k){
    auto coppia=risultato.second.front();
    gsl_matrix_set(evoluzioneAutoma,0,k,coppia.first);
    gsl_matrix_set(evoluzioneAutoma,1,k,coppia.second);
    //printf("istante: %.12lf\tstato: %lu\n",coppia.first,coppia.second);
    risultato.second.pop();
  }
  
  FILE* fileSimulazione=fopen("datiSimulazione","wb");
  fwrite(&info,sizeof(struct InfoSimulazione),1,fileSimulazione);
  gsl_matrix_fwrite(fileSimulazione,simulazione);
  gsl_matrix_fwrite(fileSimulazione,evoluzioneAutoma);
  fclose(fileSimulazione);
  
  gsl_matrix_free(risultato.first);
  gsl_matrix_free(evoluzioneAutoma);
  return 0;
}

bool Condizione1(double t, gsl_vector* y){
  double y2=gsl_vector_get(y,1);
  return y2 < -5.0;
}

bool Condizione2(double t, gsl_vector* y){
  double y1=gsl_vector_get(y,0);
  return y1 < 2.0;
}

bool Condizione3(double t, gsl_vector* y){
  double y1=gsl_vector_get(y,0);
  return y1 > 2.1;
}

void Sistema1(double t,gsl_vector* y,gsl_vector* dy){
  gsl_vector_memcpy(dy,y);
  gsl_vector_scale(dy,2.0);
}

void Sistema2(double t,gsl_vector* y,gsl_vector* dy){
  gsl_vector_memcpy(dy,y);
  gsl_vector_scale(dy,-2.0);
}

void Sistema3(double t,gsl_vector* y,gsl_vector* dy){
  gsl_vector_set_zero(dy);
}

gsl_matrix* metodoPerODE(struct InfoBaseSimulazione* info,gsl_matrix* innesco){
  double A[]={1.0};
  double B[]={0.0};
  double b_1=1.0;
  return LMM(info,A,B,b_1,innesco);
}

gsl_matrix* metodoPerInnesco(struct InfoBaseSimulazione* info,gsl_vector* x0){
  return Heun(info,x0);
}
