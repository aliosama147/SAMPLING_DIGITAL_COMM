clc
clear
close all

%% FREQUENCIES

fc = 100;
fm = fc/10;
fs = 10*fc;

%% AMPLITUDES

Amt = 1;
Act = 0.5;

%% TIME INTERVAL

t = 0 : 1/fs : 4/fm;

%% INPUT SIGNAL

mt = Amt*sin(2*pi*fm*t);

%% CARRIER SIGNAL

ct = Act*square(2*pi*fc*t)+0.5;

%% MODULATED SIGNAL(PAM)

st = mt.*ct;

%% INPUT SIGNAL FIGURE
figure
subplot(321)
plot(t,mt,'k','linewidth',2);
title('Input Signal')
xlabel('Time')
ylabel('Amplitude')

%% CARRIER SIGNAL FIGURE

subplot(322)
plot(t,ct,'b','linewidth',2);
title('Carrier Signal')
xlabel('Time')
ylabel('Amplitude')
axis([0 0.1 0 1])

%% MODULATED SIGNAL WITH NATURAL SAMPLING(PAM)

subplot(323)
plot(t,st,'r','linewidth',2);
title('Natural Sampling')
xlabel('Time')
ylabel('Amplitude')

%% MODULATED SIGNAL WITH IDEAL SAMPLING(PAM)

subplot(324)
stem(t,mt)
title('Discrete Sequence of The Signal (Ideal Sampling)')
xlabel('Time')
ylabel('Amplitude')

%% FLAT-TOP PAM CALCULATIONS

y1 = zeros(1,length(t));

for i = 2:length(t)
    if ct(i) == 1 && ct(i-1) == 0
        y1(i) = ct(i) * mt(i); % Sampling occurs during rising edge
    elseif ct(i) == 1 && ct(i-1) == 1
        y1(i) = y1(i-1); % Value remains constant while the carrier signal is 1
    else
        y1(i) = 0; % Otherwise, y is zero
    end
end

%% MODULATED SIGNAL WITH FLAT-TOP SAMPLING(PAM)

subplot(325)
plot(t, y1,'linewidth',2);
xlabel('Time');
ylabel('Amplitude');
title('Flat-Top Sampling');

%% QUANTIZATION PART


sig=mt;
partition = -1:0.3:1; % Length 11, to represent 12 intervals
codebook = -1.3:0.3:1; % Length 12, one entry for each interval
[partition2,codebook2] = lloyds(sig,codebook);
[index,quants,distor] = quantiz(sig,partition,codebook);
[index2,quant2,distor2] = quantiz(sig,partition2,codebook2);

subplot(326)
plot(t,sig,t,quant2,'-r','linewidth',2)
xlabel('Time');
ylabel('Amplitude');
title('Quantization');
legend('Original signal','Quantized signal');
axis([0 0.4 -1.2 1.2])

figure
subplot(211)
plot(t,sig,t,quants,'-r','linewidth',2)
xlabel('Time');
ylabel('Amplitude');
title('Quantization Before Optimization');
legend('Original signal','Quantized signal');
axis([0 0.4 -1.2 1.2])

subplot(212)
plot(t,sig,t,quant2,'-r','linewidth',2)
xlabel('Time');
ylabel('Amplitude');
title('Quantization After Optimization');
legend('Original signal','Quantized signal');
axis([0 0.4 -1.2 1.2])

% clc;
% p=input('Enter the probabilities:');
% n=length(p);
% symbols=[1:n];
% [dict,avglen]=huffmandict(symbols,p);
% temp=dict;
% t=dict(:,2);
% for i=1:length(temp)
%     temp{i,2}=num2str(temp{i,2});
% end
% disp('The huffman code dict:');
% disp(temp)
% fprintf('Enter the symbols between 1 to %d in[]',n);
% sym=input(':')
% encod=huffmanenco(sym,dict);
% disp('The encoded output:');
% disp(encod);
% bits=input('Enter the bit stream in[];');
% decod=huffmandeco(bits,dict);
% disp('The symbols are:');
% disp(decod); 
% 
% H=0;
% Z=0;
% for(k=1:n)
%     H=H+(p(k)*log2(1/p(k)));
%         
%  end
% fprintf(1,'Entropy is %f bits',H);
% N=H/avglen;
% fprintf('\n Efficiency is:%f',N);
% for(r=1:n)
%    l(r)=length(t{r});
% end
% m=max(l)
% s=min(l)
% v=m-s;
% fprintf('the variance is:%d',v);
