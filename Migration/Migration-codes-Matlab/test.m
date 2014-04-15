dd=zeros(256,64);
dd(50,32)=1;
dd(100,32)=-1;
dd(200,32)=1;
data=conv2(dd,hamming(5));
data=data(1:256,:);

 M=GAZDAG_adj(data(1:256,:),.004,10*(1:64),5,5,1000,1000*ones(1,200),5,60,4,4);
 figure(1)
 subplot(1,2,1)
 wigb(data,.1,10*(1:64),(0:255)*.004)
 xlabel('CMP location (m)')
 ylabel('Time (s)')
 title('Spike data for test')
 subplot(1,2,2)
 wigb(M,1,10*(1:64),5*(1:200))
 xlabel('CMP location (m)')
 ylabel('Depth (m)')
 title('Impulse responce test')
 
 [d,h,t] = linear_events(.004,30,1,(0:100)*10,[.2 .3 .5],[1/1500,1/2500,-1/2000],[1,-1,1],1000,1);
 
 M=GAZDAG_adj(d,.004,h,5,5,1000,1000*ones(1,200),5,60,4,4);
 
 figure(2)
 subplot(1,2,1)
 wigb(d,1,(0:100)*10,(0:length(d(:,1))-1)*.004)
xlabel('CMP location (m)')
 ylabel('Time (s)')
  title('Linear data')
 subplot(1,2,2)
 wigb(M(:,1:100),1,(0:99)*10,5*(1:200))
 xlabel('CMP location (m)')
 ylabel('Depth (m)')
  title('migrated section of Linear data')
  
 D=GAZDAG_forward(M,10,(0:199)*5,0.004,1,1000*ones(200,1),4,4);
   figure(3)
 subplot(1,2,1)
 wigb(d,1,(0:100)*10,(0:length(d(:,1))-1)*.004)
xlabel('CMP location (m)')
 ylabel('Time (s)')
  title('Linear data')
 subplot(1,2,2)
 wigb(D(:,1:100),1,(0:99)*10,(0:length(d(:,1))-1)*.004)
 xlabel('CMP location (m)')
 ylabel('Depth (m)')
  title('GAZDAG de- migrated data')