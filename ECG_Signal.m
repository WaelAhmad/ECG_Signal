ecg = load('DataN.txt');
fs = 256;
time_period = 1/fs;
x = (1/(8*time_period));
nyquist_rate = fs/2;
w0 = 50/nyquist_rate;
bw = w0/10;


[b,a] = iirnotch(w0,bw);
notchdata= filter(b,a,ecg);
nsamples=size(notchdata,1);
time=zeros(nsamples,1);
 for i = 1:1:nsamples+1
     time(i)=i/256;
 end
figure(1);
plot(1:2000,ecg(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal before filtering');

figure(2);
plot((1:2000),notchdata(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal after Notch Filter');

low = 0.1 / nyquist_rate;
high = 45 / nyquist_rate;
wn = [low high];
[b,a] = butter(2,wn,'bandpass');
banddata= filter(b,a,notchdata);
nsamples=size(banddata,1);
figure(3);
plot((1:2000),notchdata(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal after BPF');

figure(4), subplot 211
plot((1:2000),notchdata(1:2000))
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal before noise filtering');
legend('BEFORE FILTERING');
subplot 212
plot((1:2000),notchdata(1:2000),'r')
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal after noise filtering');
legend('AFTER FILTERING');

derivativedata = zeros(nsamples,1);
for n = 3:1:nsamples-3
        derivativedata(n) = x * (-banddata(n-2) - 2*banddata(n-1) + 2*banddata(n+1) + banddata(n+2));
end
figure(5);
plot((1:2000),derivativedata(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal after differentiation');

squaredata = derivativedata.*derivativedata;
figure(6);
plot((1:2000),squaredata(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal after squaring');

windowSize = 25;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
avgmovingdata = filter(b,a,squaredata);
figure(7);
plot((1:2000),avgmovingdata(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal after average moving window');

threshold = 0.6 * max(avgmovingdata);
N1=length(avgmovingdata);
j=[];
v=[];
m=1;
n = windowSize;
n = 25;
i=1;
while n<=2000
    y = max(avgmovingdata(m:n));
    index=find(avgmovingdata==y);
    
    if(y>threshold) 
        j(i)=index;
        v(i)=y;
        i=i+1;
        
    end
    n=n+windowSize;
    m=m+windowSize-1;
end
size(j)
figure(8);
plot(avgmovingdata(1:2000));
title('ECG signal the detected R waves for N=25');
hold on
plot(j,v,'*')
title('ECG signal the detected R waves for N=25');
xlabel('Samples');
ylabel('Amplitude');



Nsample=size(ecg,1);
derivativedata = zeros(Nsample,1);
for n = 3:1:Nsample-3
        derivativedata(n) = x * (-ecg(n-2) - 2*ecg(n-1) + 2*ecg(n+1) + ecg(n+2));
end
% figure(9);
% plot((1:2000),derivativedata(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal after differentiation  without noise filtering');

squaredata = derivativedata.*derivativedata;
% figure(10);
% plot((1:2000),squaredata(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal the detected R waves for N=25');

windowSize = 25;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
avgmovingwithnoise = filter(b,a,squaredata);
figure(11);
plot((1:2000),avgmovingwithnoise(1:2000));
xlabel('Samples');
ylabel('Amplitude');
title('ECG signal after average moving window without noise filtering');

findpeaks(avgmovingwithnoise(1:2000),'MinPeakHeight',100,'MinPeakDistance',200);
axis([0 2000 800 2000]);
title('ECG signal without noise filtering the detected R waves for N=25');


% threshold1 = 0.6 * max(avgmovingwithnoise);
% N=length(avgmovingwithnoise);
% u=[];
% s=[];
% n = windowSize;
% n = 25;
% f=1;
% e=1;
% 
% 
% while n<=2000
%     k = max(avgmovingwithnoise(f:n));
%     index=find(avgmovingwithnoise==k);
%     
%     if(k>threshold1) 
%         u(e)=index;
%         s(e)=k;
%         e=e+1;
%         
%     end
%     n=n+windowSize;
%     f=f+windowSize-1;
% end
% size(u)
% figure(12);
% plot(avgmovingwithnoise(1:2000));
% hold on
% plot(u,s,'*r')
% title('ECG signal without noise filtering the detected R waves for N=25');
% xlabel('Samples');
% ylabel('Amplitude');

% p=max(avgmovingwithnoise);
% threshold1=0.6*max(avgmovingwithnoise);
% N1=length(avgmovingwithnoise);
% u=zeros(N1,1);
% for s=2:N1-1
%     if((avgmovingwithnoise(s)>avgmovingwithnoise(s+1))&&(avgmovingwithnoise(s)>avgmovingwithnoise(s-1))&&(avgmovingwithnoise(s)>threshold1)) 
%         u(s,1)=s;
%     end
% end
% figure(12);
% plot(avgmovingwithnoise(1:2000));
% hold on;
% plot(u(1:2000),p,'*r');



