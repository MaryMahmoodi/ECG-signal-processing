function [ locs, feature_vector ] = ECG_R_QRS_detection( ecg,fs , show)

%output
% feature_vector=[QRS_max QRS_ave QRS_min Heartrate];
%QRS_max: maximum duration of QRS in the ECG signal
% QRS_min: minimum duration of QRS in the ECG signal 
% locs( sample for R location), QRS_max (ms), QRS_ave(ms), QRS_min(ms), Heartrate(beats/s) 

% Example: 
%load ecg_hfn1;fs=1000; 
%[ locs, feature_vector ] = ECG_R_QRS_detection( ecg_hfn1,1000 , 1);



% if show; it shows the figures
if nargin<3;show=0;end
slen = length(ecg);
t=[1:slen]/fs;

if show
   figure(1); clf
subplot(4,1,1);plot(t,ecg); title('original ECG signal')
axis tight;
end
ecg=ecg-mean(ecg);



% % % % % d1= fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2,C',40,0.01,0.16,10,12,fs);%50,0.01,0.16,45,50.5,fs);%36
% % % % % Hd1=design(d1,'equiripple');
% % % % % ecg1=filter(Hd1,ecg);

% if fs>500
    [a1,b1]=butter(3,0.01/fs/2,'high');
% % % [a11,b11]=butter(3,30/fs/2,'low');
 ecg1=filtfilt(a1,b1,ecg);%filter
% % %  ecg1=filtfilt(a11,b11,ecg1);
% else
%     ecg1=ecg;
% end
ecg=ecg-mean(ecg);

for i=1:length(ecg1)
    if isnan(ecg1(i));
        ecg1(i)=0;
    end
end

if (ecg1)==zeros(size(ecg1))
    ecg=ecg;
else
ecg=ecg1;
end
if show
subplot(4,1,2);plot(t,ecg); title('Highpass filtered ECG signal')
axis tight;
end
%  [pks,locs] = findpeaks(ecg_m,'MINPEAKHEIGHT',0.0035);%0.0035 try and error it is used for s1
for i=1:length(ecg)
    if isnan(ecg(i))
        ecg(i)=0;
    end
end
% % % % % % % % % %   [ecg, cost] = tvd_mm(ecg, 4, 5);


%% spike detection
threshold=abs( 1*1/length(ecg)*sum((ecg))+3*std((ecg)) ); %2 is good
spikes=zeros(size(ecg));

    for i=1:length(ecg)
            if ecg(i)>threshold
                spikes(i)=ecg(i);
            else
                                spikes(i)=0;
            end
    end

org_spikes=spikes;
if show
ax(3)=subplot(413);plot(t,org_spikes);axis tight;title('R peaks');
end

binary =abs(org_spikes)>0;
% find consequtive detected spike samples
E = binary(2:end)-binary(1:end-1);
sise = size(binary);

begins = find(E==1)+1;

if binary(1) == 1
    if sise(1) > 1
        begins = [1; begins];
    elseif sise(2) > 1
        begins = [1 begins];
    else
        error('The input signal is not one dimensional')
    end
elseif numel(begins) == 0 && binary(1) == 0
    begins = NaN;
end

ends = find(E==-1);
if binary(end) == 1
    if sise(1) > 1
        ends = [ends; length(binary)];
    elseif sise(2) > 1
        ends = [ends length(binary)];
    else
        error('The input signal is not one dimensional')
    end
elseif numel(ends) == 0 && binary(end) == 0
    ends = NaN;
end
ecg1=abs(ecg);
locs=[];QRS=[];interval=[];R=[];S=[];Q=[];
if ~isnan(begins)
for i=1:length(begins)
    R=[];S=[];Q=[];
    R=find(ecg==max(ecg(begins(i) : ends(i))));
    if R(1)+0.05*fs > length(ecg)
            S=find(ecg1==min( ecg1( R(1):R(1)+0.0*fs ) ) )  ;
    else
    S=find(ecg1==min( ecg1( R(1):R(1)+0.05*fs ) ) )  ;
    end
    if R(1)-0.05*fs > 0
    Q=find( ecg1==min(ecg1(R(1)-0.05*fs:R(1)) )  )  ;
    else
            Q=find( ecg1==min(ecg1(R(1)-0.0*fs:R(1)) )  )  ;
    end
    
    locs(1,i)=R(1);%round(mean([begins(i) ends(i) ]));% sample for R location
  QRS(1,i)= abs((S(1)- Q(end)))/fs *1000;%abs(begins(i) - ends(i))/fs*1000;%abs((S(1)- Q(end)))/fs *1000;
  
  if i>1
      interval(1,i-1)=locs(i)-locs(i-1);
  end
  end
end
% locs( sample for R location), QRS_max (ms), QRS_ave(ms), QRS_min(ms), Heartrate(beats/s) 
if ~ isempty (QRS)
QRS_max=max(QRS);
QRS_ave=mean(QRS);
QRS_min=min(QRS);
if QRS_min==0
    QRS_min=1;
end
else
 QRS_max=1;
 QRS_ave=1;
 QRS_min=1;
end

if isempty (interval);
    Heartrate=0;
else
Heartrate=60/(max(interval)/fs);
end
feature_vector=[QRS_max QRS_ave QRS_min Heartrate];
% locs( sample for R location), QRS_max (ms), QRS_ave(ms), QRS_min(ms), Heartrate(beats/s) 


if show
ax(4)=subplot(414);hold on
plot(t,ecg);hold on;
plot(t(locs),ecg(locs),'k^','markerfacecolor',[1 0 0]);title('OUTPUT R-peaks');
 axis tight;
end


end

