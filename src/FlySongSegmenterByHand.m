function FlySongSegmenterByHand(data, Fs, varargin)
%function FlySongSegmenterByHand(data,Fs, pulseTimes)
%
%data is a matrix containing the data (TxNChannels)
%Fs is the sampling rate (Hz)
%pulseTimes (optional) vector of predefined pulseTimes (in samples)
%
%popupmenu on the left selects the channel number
%pan with < and >, zoom with ^ and v
%`p` activates pulse annotation mode: 
%  add multiple pulses by clicking on each and then pressing
%  the return key (PULSE_MODE=1)
%  add a pulse by clicking on the peak and trough (PULSE_MODE=2)
%add a sine song by clicking on the beginning and end
%delete removes just those currently displayed
%save creates a .mat file with _byhand appended to filename
%  if workspace data passed in, file is workspace_byhand.mat


global RAW IDXP IDXS XPAN XZOOM YPAN YZOOM UNITS TYPE
global CHANNEL DATA PULSE SINE PULSE_MODE NFFT
global FS NARGIN H BISPECTRUM FTEST PARAMS NW K TEXT
global SHOWPULSESFROMALLCHANNELS SHOWTRACESFROMALLCHANNELS


UNITS{1}=1;
UNITS{2}='s';

PULSE_MODE=1;
SHOWPULSESFROMALLCHANNELS = true; % display pulses for all channels in trace or not
SHOWTRACESFROMALLCHANNELS = true; % display all traces in one plot - for multichannel recordings
IDXP=[];
IDXS=[];
XPAN=0;  % sec
XZOOM=0.1;  % sec
CHANNEL=1;
NARGIN=nargin;
NFFT=2^9;  % tic
BISPECTRUM=0;
FTEST=0;
PARAMS=[];
NW=9;  K=17;


TYPE='';
PULSE=[];
SINE=[];
DATA=data;
NCHAN=size(DATA,2);
FS=Fs;
RAW=DATA(:,CHANNEL);

if nargin==3
   PULSE=varargin{1};
   if size(PULSE,1)==1
      PULSE = PULSE';
   end
   if size(PULSE,2)==1
      PULSE = [ones(size(PULSE)), PULSE];
   end
end

YPAN=0;  % Hz
YZOOM=FS/2;  % Hz

if(XPAN>length(RAW)/FS) 
   XPAN=0; 
end
if((XPAN+XZOOM)>(length(RAW)/FS))
   XZOOM=(length(RAW)/FS)-XPAN; 
end

figure;
tmp=get(gcf,'position');
set(gcf,'position',[0 0 2*tmp(3) 1.5*tmp(4)]);
set(gcf,'menubar','none','ResizeFcn',@resizeFcn,'WindowKeyPressFcn',@windowkeypressFcn);

H=uipanel();
AX = axes('Units','normal', 'Position', [0 0 1 1], 'Parent', H);

uicontrol('parent',H,'style','popupmenu','value',CHANNEL,...
   'string',1:NCHAN, ...
   'callback', @changechannel_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','-','tooltipstring','pan X left 10x', ...
   'callback', @panleft10x_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','<','tooltipstring','pan X left; right click is zoom Y out', ...
   'callback', @panleft_callback, ...
   'buttondownfcn', @yzoomout_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','^','tooltipstring','zoom X in;  right click is pan Y up', ...
   'callback', @xzoomin_callback,...
   'buttondownfcn', @panup_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','v','tooltipstring','zoom X out;  right click is pan Y down', ...
   'callback', @xzoomout_callback,...
   'buttondownfcn', @pandown_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','>','tooltipstring','pan X right;  right click is zoom Y in', ...
   'callback', @panright_callback,...
   'buttondownfcn', @yzoomin_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','+','tooltipstring','pan X right 10x', ...
   'callback', @panright10x_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(f)req','tooltipstring','increase frequency resolution', ...
   'callback', @nfftup_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(t)ime','tooltipstring','increase temporal resolution', ...
   'callback', @nfftdown_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(b)ispectrum','tooltipstring','toggle between spectrum and bispectrum', ...
   'callback', @bispectrum_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(m)ulti-taper F-test','tooltipstring','toggle between spectrum and multi-taper spectrum with F-test in green;  use N/n and K/k to change bandwidth and # tapers respectively', ...
   'callback', @ftest_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(7) brown & puckette','tooltipstring','track chosen harmonic with high resolution', ...
   'callback', @brown_puckette_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(p)ulse','tooltipstring','add pulse song', ...
   'callback', @addpulse_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(s)ine','tooltipstring','add sine song', ...
   'callback', @addsine_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(d)elete','tooltipstring','delete displayed pulse and sine song', ...
   'callback', @delete_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','(l)isten','tooltipstring','listen to displayed recording', ...
   'callback', @listen_callback);
uicontrol('parent',H,'style','pushbutton',...
   'string','save','tooltipstring','save segmentation to disk', ...
   'callback', @save_callback);
TEXT=uicontrol('parent',H,'style','text',...
   'horizontalalignment','left',...
   'tooltip','time and frequency resolution of FFT (and time-bandwidth product and # of tapers for F-test)');

update;


function resizeFcn(src,evt)
global H

tmp=get(gcf,'position');
recording=get(H,'children');

for i = 1:length(recording)
  if strfind(class(recording(i)), 'Axes')
     continue
  end
  if(strcmp(get(recording(i),'style'),'popupmenu'))
    set(recording(i),'position',[10,tmp(4)-30,70,20]);
  elseif(strcmp(get(recording(i),'style'),'text'))
    set(recording(i),'position',[10,10,350,15]);
  else
    switch(get(recording(i),'string'))
      case('-')
        set(recording(i),'position',[80,tmp(4)-30,20,20]);
      case('<')
        set(recording(i),'position',[100,tmp(4)-30,20,20]);
      case('^')
        set(recording(i),'position',[120,tmp(4)-30,20,20]);
      case('v')
        set(recording(i),'position',[140,tmp(4)-30,20,20]);
      case('>')
        set(recording(i),'position',[160,tmp(4)-30,20,20]);
      case('+')
        set(recording(i),'position',[180,tmp(4)-30,20,20]);
      case('(f)req')
        set(recording(i),'position',[200,tmp(4)-30,40,20]);
      case('(t)ime')
        set(recording(i),'position',[240,tmp(4)-30,40,20]);
      case('(b)ispectrum')
        set(recording(i),'position',[280,tmp(4)-30,80,20]);
      case('(m)ulti-taper F-test')
        set(recording(i),'position',[360,tmp(4)-30,120,20]);
      case('(7) brown & puckette')
        set(recording(i),'position',[480,tmp(4)-30,120,20]);
      case('(p)ulse')
        set(recording(i),'position',[600,tmp(4)-30,50,20]);
      case('(s)ine')
        set(recording(i),'position',[650,tmp(4)-30,50,20]);
      case('(d)elete')
        set(recording(i),'position',[700,tmp(4)-30,50,20]);
      case('(l)isten')
        set(recording(i),'position',[750,tmp(4)-30,50,20]);
      case('save')
        set(recording(i),'position',[800,tmp(4)-30,40,20]);
    end
  end
end


function windowkeypressFcn(src,evt)
global NW K NFFT

switch(evt.Key)
  case('hyphen')
    panleft10x_callback;
  case('leftarrow')
    if(strcmp(evt.Modifier,'shift'))
      yzoomout_callback;
    else
      panleft_callback;
    end
  case('uparrow')
    if(strcmp(evt.Modifier,'shift'))
      panup_callback;
    else
      xzoomin_callback;
    end
  case('downarrow')
    if(strcmp(evt.Modifier,'shift'))
      pandown_callback;
    else
      xzoomout_callback;
    end
  case('rightarrow')
    if(strcmp(evt.Modifier,'shift'))
      yzoomin_callback;
    else
      panright_callback;
    end
  case('equal')
    panright10x_callback;
  case('f')
    nfftup_callback;
  case('t')
    nfftdown_callback;
  case('b')
    bispectrum_callback;
  case('m')
    ftest_callback;
  case('7')
    brown_puckette_callback;
  case('p')
    addpulse_callback;
  case('s')
    addsine_callback;
  case('d')
    delete_callback;
  case('l')
    listen_callback;
  case('n')
    if(strcmp(evt.Modifier,'shift'))
      if((NW+1)<(NFFT/2))
        NW=NW+1;
      end
    else
      NW=max(ceil((K+1)/2),NW-1);
    end
    update;
  case('k')
    if(strcmp(evt.Modifier,'shift'))
      K=min(2*NW-1,K+1);
    else
      K=K-1;
    end
    update;
end



function changechannel_callback(hObject,eventdata)
global DATA RAW CHANNEL

recording=CHANNEL;
CHANNEL=get(hObject,'value');
if(recording==CHANNEL)  beep;  return;  end
RAW=DATA(:,CHANNEL);
update;



function panleft_callback(hObject,eventdata)
global XPAN XZOOM;

recording=XPAN;
XPAN=max(0,XPAN-XZOOM/2);
if(recording~=XPAN)  update;  else  beep;  end



function panright_callback(hObject,eventdata)
global XPAN XZOOM RAW FS;

recording=XPAN;
XPAN=min(length(RAW)/FS-XZOOM,XPAN+XZOOM/2);
if(recording~=XPAN)  update;  else  beep;  end



function panleft10x_callback(hObject,eventdata)
global XPAN XZOOM;

recording=XPAN;
XPAN=max(0,XPAN-5*XZOOM);
if(recording~=XPAN)  update;  else  beep;  end



function panright10x_callback(hObject,eventdata)
global XPAN XZOOM RAW FS;

recording=XPAN;
XPAN=min(length(RAW)/FS-XZOOM,XPAN+5*XZOOM);
if(recording~=XPAN)  update;  else  beep;  end



function panup_callback(hObject,eventdata)
global YPAN YZOOM RAW FS;

recording=YPAN;
YPAN=min(FS/2-YZOOM,YPAN+YZOOM/2);
if(recording~=YPAN)  update;  else  beep;  end



function pandown_callback(hObject,eventdata)
global YPAN YZOOM RAW FS;

recording=YPAN;
YPAN=max(0,YPAN-YZOOM/2);
if(recording~=YPAN)  update;  else  beep;  end



function xzoomin_callback(hObject,eventdata)
global XPAN XZOOM FS NFFT; 

%if(XZOOM*FS<10*NFFT)  return;  end;
recording=XPAN;  bar=XZOOM;
XZOOM=XZOOM/2;
XPAN=XPAN+XZOOM/2;
if((recording~=XPAN)||(bar~=XZOOM))  update;  else  beep;  end



function yzoomin_callback(hObject,eventdata)
global YPAN YZOOM FS NFFT;

if(YZOOM<10*FS/NFFT)  return;  end;
recording=YPAN;  bar=YZOOM;
YZOOM=YZOOM/2;
YPAN=YPAN+YZOOM/2;
if((recording~=YPAN)||(bar~=YZOOM))  update;  else  beep;  end



function xzoomout_callback(hObject,eventdata)
global XPAN XZOOM RAW FS;

recording=XPAN;  bar=XZOOM;
XPAN=max(0,XPAN-XZOOM/2);
XZOOM=XZOOM*2;
if((XPAN+XZOOM)>(length(RAW)/FS))
  XZOOM=(length(RAW)/FS)-XPAN;
end
if((recording~=XPAN)||(bar~=XZOOM))  update;  else  beep;  end



function yzoomout_callback(hObject,eventdata)
global YPAN YZOOM RAW FS;

recording=YPAN;  bar=YZOOM;
YPAN=max(0,YPAN-YZOOM/2);
YZOOM=YZOOM*2;
if((YPAN+YZOOM)>(FS/2))
  YZOOM=FS/2-YPAN;
end
if((recording~=YPAN)||(bar~=YZOOM))  update;  else  beep;  end



function addpulse_callback(hObject,eventdata)
global CHANNEL PULSE PULSE_MODE;

if(PULSE_MODE==1)
  tmp=ginput;
  tmp2=size(tmp,1);
  PULSE(end+1:end+tmp2,:)=[repmat(CHANNEL,tmp2,1) tmp(:,1)];
else
  tmp=ginput(2);
  PULSE(end+1,:)=[CHANNEL tmp(:,1)'];
end
update;



function addsine_callback(hObject,eventdata)
global CHANNEL SINE;

tmp=ginput(2);
SINE(end+1,:)=[CHANNEL tmp(:,1)'];
update;



function delete_callback(hObject,eventdata)
global PULSE SINE IDXP IDXS;

tmpp=setdiff(1:size(PULSE,1),IDXP);
tmps=setdiff(1:size(SINE,1),IDXS);
recording='yes';
bar=length(IDXP)+length(IDXS);
if(bar>10)
   recording=questdlg(['are you sure you want to delete these ' num2str(bar) ' items?'],...
       '','yes','no','no');
end
if(strcmp(recording,'yes'))
  PULSE=PULSE(tmpp,:);
  SINE=SINE(tmps,:);
end
update;



function listen_callback(hObject,eventdata)
global RAW XPAN XZOOM FS;

sound(RAW((1+ceil(XPAN*FS)):floor((XPAN+XZOOM)*FS)),min(48000,FS));



function save_callback(hObject,eventdata)
global PULSE SINE NARGIN

if(NARGIN==1)
  save(['_byhand.mat'],'PULSE','SINE');
else
  save(['workspace_byhand.mat'],'PULSE','SINE');
end



function nfftup_callback(hObject,eventdata)
global NFFT NW K FTEST FS XZOOM

if(XZOOM*FS<10*NFFT)  return;  end;
NFFT=NFFT*2;
update;



function nfftdown_callback(hObject,eventdata)
global NFFT NW K FTEST FS YZOOM

if(YZOOM<10*FS/NFFT)  return;  end;
recording=NFFT;
NFFT=max(8,NFFT/2);
if(FTEST)
  while(NW>=(NFFT/2))
    NW=NW-1;
  end
  if(K>(2*NW-1))
    K=floor(2*NW-1);
  end
end
if(recording~=NFFT)  update;  else  beep;  end



function bispectrum_callback(hObject,eventdata)
global BISPECTRUM FTEST

if(FTEST)  return;  end

BISPECTRUM=~BISPECTRUM;
update;



function ftest_callback(hObject,eventdata)
global FTEST BISPECTRUM NW K NFFT FS

if(BISPECTRUM)  return;  end

FTEST=~FTEST;
if(FTEST)
  while(NW>(NFFT/2))
    NW=NW-1;
  end
  if(K>(2*NW-1))
    K=floor(2*NW-1);
  end
end
update;



function brown_puckette_callback(hObject,eventdata)
global FS RAW XPAN XZOOM NFFT

[t,f]=ginput(1);

recording=RAW((round(t*FS)-NFFT/2):floor((XPAN+XZOOM)*FS));
[Hr tr]=brown_puckette(recording',NFFT,FS,1,f);

T=flipud(RAW((1+ceil(XPAN*FS)):(round(t*FS)+NFFT/2-1)));
[Hl tl]=brown_puckette(T',NFFT,FS,1,f);

H=[fliplr(Hl) Hr];
t=[t-fliplr(tl)+tl(1) t+tr-tr(1)];

plot3(t,H,ones(size(t)),'r.-');


function update()
global RAW DATA IDXP IDXS XPAN XZOOM YPAN YZOOM UNITS
global CHANNEL PULSE SINE PULSE_MODE NFFT
global FS BISPECTRUM FTEST PARAMS NW K TEXT
global SHOWPULSESFROMALLCHANNELS SHOWTRACESFROMALLCHANNELS

if(FTEST)
  set(TEXT,'string',['dT=' num2str(NFFT/FS) ' s, dF=' num2str(FS/NFFT) ' Hz, NW=' num2str(NW) ' (' num2str(round(NW/NFFT*FS)) ' Hz), K=' num2str(K)]);
else
  set(TEXT,'string',['dT=' num2str(NFFT/FS) ' s, dF=' num2str(FS/NFFT) ' Hz']);
end

T=(1+ceil(XPAN*FS)):floor((XPAN+XZOOM)*FS);
recording=RAW(T);

subplot(2,1,1);  cla;  hold on;
if(XZOOM*FS>10*NFFT)
if(length(recording)>NFFT)
  if(FTEST)
    if(isempty(PARAMS)||(PARAMS.NW~=NW)||(PARAMS.K~=K)||(PARAMS.NFFT~=NFFT))
      PARAMS=[];
      PARAMS.NW=NW;
      PARAMS.K=K;
      PARAMS.NFFT=NFFT;
      [PARAMS.tapers,~]=dpsschk([NW K],NFFT,FS);
      PARAMS.Fs=FS;
      PARAMS.pad=0;
      PARAMS.fpass=[0 FS/2];
    end
    kk=ceil((length(recording)-NFFT+1)/(NFFT/2));
    t=(0:(kk-1))*NFFT/2/FS+NFFT/2/FS;

    tmp=2^nextpow2(NFFT)/2+1;
    p=zeros(tmp,kk);
    F=zeros(tmp,kk);
    f=zeros(tmp,kk);
    sig=zeros(kk);
    sd=zeros(tmp,kk);
    t=(0:(kk-1))*NFFT/2/FS+NFFT/2/FS;
    parfor k=1:kk
      i=1+(k-1)*NFFT/2;
      [F(:,k),p(:,k),f(:,k),sig(k),sd(:,k)] = ftestc(recording(i:(i+NFFT-1)),PARAMS,0.01/NFFT,'n');
    end
    f=f(:,1);  sig=sig(1);  sd=sd(:,1);

  else
    [s,f,t,p]=spectrogram(recording',NFFT,[],[],FS,'yaxis');
    if(BISPECTRUM)
      [p,f]=my_polyspectrum(s',f,[1 2 -3]);
      p=p';
    end
  end
  fidx=find(f>=YPAN & f<=(YPAN+YZOOM));
  tmp=log10(abs(p(fidx,:)));
  tmp4=reshape(tmp,1,prod(size(tmp)));
  tmp2=prctile(tmp4,1);
  tmp3=prctile(tmp4,99);
  idx=find(tmp<tmp2);  
  tmp(idx)=tmp2;
  idx=find(tmp>tmp3);  
  tmp(idx)=tmp3;
  surf((t+T(1)./FS).*UNITS{1},f(fidx)-f(2)/2,tmp,'EdgeColor','none');
  if(BISPECTRUM)
    colormap(gray);
  else
    colormap(flipud(gray));
    if(FTEST)
      [i,j]=find(F>repmat(sig,size(F)));
      tmp=t(2)-t(1);  tmp2=f(2)-f(1);
      for(k=1:length(i))
        plot((t(j(k))+[0 tmp tmp 0 0]+T(1)/FS),...
            f(i(k))-f(2)/2+[0 0 tmp2 tmp2 0],'g-');
      end
    end
  end
  axis tight;
  ylabel('frequency (Hz)');
end
end

subplot(2,1,2);  cla;  hold on;

IDXP=[];
if(~isempty(PULSE))
  if(PULSE_MODE==1)
     if SHOWPULSESFROMALLCHANNELS% plot pulses in all channels
         IDXP=find((PULSE(:,2)>=XPAN) & (PULSE(:,2)<=(XPAN+XZOOM)));
     else
         IDXP=find((PULSE(:,1)==CHANNEL) & ...
                 (((PULSE(:,2)>=XPAN) & (PULSE(:,2)<=(XPAN+XZOOM)))));
     end
  else
    IDXP=find((PULSE(:,1)==CHANNEL) & ...
           (((PULSE(:,2)>=XPAN) & (PULSE(:,2)<=(XPAN+XZOOM))) | ...
            ((PULSE(:,3)>=XPAN) & (PULSE(:,3)<=(XPAN+XZOOM)))));
  end
end

IDXS=[];
if(~isempty(SINE))
  IDXS=find((SINE(:,1)==CHANNEL) & ...
           (((SINE(:,2)>=XPAN) & (SINE(:,2)<=(XPAN+XZOOM))) | ...
            ((SINE(:,3)>=XPAN) & (SINE(:,3)<=(XPAN+XZOOM)))));
end

if SHOWTRACESFROMALLCHANNELS
   axisLimits = [min(min(DATA(T,:))) max(max(DATA(T,:)))  max(max(DATA(T,:)))  min(min(DATA(T,:)))];
else
   axisLimits = [min(recording)  max(recording)  max(recording)  min(recording)];
end

for(i=1:length(IDXP))
  if(PULSE_MODE==1)
    plot([PULSE(IDXP(i),2) PULSE(IDXP(i),2)],axisLimits(1:2),...
          'b-','linewidth',3);
  else
    patch([PULSE(IDXP(i),2) PULSE(IDXP(i),2) PULSE(IDXP(i),3) PULSE(IDXP(i),3)],...
            axisLimits, 'b','EdgeColor','b');
  end
end
for(i=1:length(IDXS))
  patch([SINE(IDXS(i),2) SINE(IDXS(i),2) SINE(IDXS(i),3) SINE(IDXS(i),3)],...
     axisLimits, 'r','EdgeColor','r');
end

if SHOWTRACESFROMALLCHANNELS
   plot(T./FS.*UNITS{1},DATA(T,:),'Color', [.6 .6 .6]);
end
plot(T./FS.*UNITS{1},recording,'k');

axis tight;
v=axis;
axis([T(1)./FS.*UNITS{1} T(end)./FS.*UNITS{1} v(3) v(4)]);
xlabel(['time (' UNITS{2} ')']);

v=axis;
subplot(2,1,1);
vv=axis;
axis([v(1) v(2) YPAN YPAN+YZOOM]);