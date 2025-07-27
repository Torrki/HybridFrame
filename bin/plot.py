import struct
import numpy as np
import matplotlib.pyplot as plt

with open("datiSimulazione","rb") as fileDati:
    dati=fileDati.read()
    dimensioneInfo=int(dati[0])
    print(dimensioneInfo)
    pattern='NNNNdd'
    (dimensioneInfo,n,N,NAutoma,t0,T)=struct.unpack(pattern+'x'*(dimensioneInfo-struct.calcsize(pattern)),dati[0:dimensioneInfo])
    print(dimensioneInfo,n,N,t0,T)
    dimensioneDati=len(dati)-dimensioneInfo
    
    O_sim=np.zeros((n,N))
    O_automa=np.zeros((2,NAutoma+1));
    cursore=dimensioneInfo
    for k in range(0,n):
        O_sim[k]=np.array( list(struct.unpack('d'*N,dati[cursore:cursore+(N*8)])) )
        cursore = cursore+(N*8)
        
    O_automa[0]=np.array( list(struct.unpack('d'*NAutoma,dati[cursore:cursore+(NAutoma*8)])) + [0.0] )
    cursore = cursore+(NAutoma*8)
    O_automa[1]=np.array( list(struct.unpack('d'*NAutoma,dati[cursore:cursore+(NAutoma*8)])) + [0.0])
    O_automa[0,NAutoma]=t0+T
    O_automa[1,NAutoma]=O_automa[1,NAutoma-1]
        
tspan=np.linspace(t0,t0+T,N)
plt.plot(tspan,O_sim.T)
#plt.step(O_automa[0],O_automa[1],where='post')
plt.grid(True)
plt.show()
