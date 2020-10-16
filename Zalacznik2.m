%PROGRAM DO WERYFIKACJI DZIA£ANIA ALGORYTMU DO TRÓJWYMIAROWEJ TRAJEKTORII ZNACZNIKA NA CIELE
%Aleksandra Michalska

clear
close all
clc

%%
zasieg = 9;                     %zasiêg znacznika pomiêdzy klatkami
znana_dlugosc = 508;

%%
%wczytywanie danych
I = vision.VideoFileReader('weryfikacja_front.avi');
Im1 = VideoReader('weryfikacja_front.avi');
videoPlayer = vision.VideoPlayer('Position',[100,100,720,480]);
nframes1 = Im1.NumberOfFrames;

I2 = vision.VideoFileReader('weryfikacja_bok.avi');
Im2 = VideoReader('weryfikacja_bok.avi');
videoPlayer2 = vision.VideoPlayer('Position',[100,100,720,480]);
nframes2 = Im2.NumberOfFrames;

%%
%macierze ze wspó³rzêdnymi œrodków markera w kolejnych klatkach
srodki1 = zeros(nframes1,2);
srodki2 = zeros(nframes2,2);

%zaznaczenie punktu bazowego na pierwszej klatce obu nagrañ
base_frame1 = I();   
figure                        
imshow(base_frame1);
title({'Pierwsza klatka nagrania I.';'1. Zaznacz obszar wokó³ znacznika.';'2. Zaznacz œrodek markera.'});
objectRegion=round(getPosition(imrect));
[base_xp1,base_yp1]=ginput(1);
srodki1(1,1)=base_xp1;
srodki1(1,2)=base_yp1;

base_frame2 = I2(); 
figure
imshow(base_frame2);
title({'Pierwsza klatka nagrania II.';'1. Zaznacz obszar wokó³ znacznika.';'2. Zaznacz œrodek markera.'});
objectRegion2=round(getPosition(imrect));
[base_xp2,base_yp2]=ginput(1);
srodki2(1,1)=base_xp2;
srodki2(1,2)=base_yp2;

close all

otoczenie1 = objectRegion(1,3)/2;
otoczenie2 = objectRegion2(1,3)/2;

figure
ref_block1=base_frame1(round(srodki1(1,2)-otoczenie1):round(srodki1(1,2)+otoczenie1),round(srodki1(1,1)-otoczenie1):round(srodki1(1,1)+otoczenie1),:);
subplot(1,2,1)
imshow(ref_block1)
title('Otoczenie œrodka znacznika nagrania I, ograniczonego przez parametr "otoczenie1"');
ref_block2=base_frame2(round(srodki2(1,2)-otoczenie2):round(srodki2(1,2)+otoczenie2),round(srodki2(1,1)-otoczenie2):round(srodki2(1,1)+otoczenie2),:);
subplot(1,2,2)
imshow(ref_block2)
title('Otoczenie œrodka znacznika nagrania II, ograniczonego przez parametr "otoczenie2"');

%%

%Pêtla wyznaczaj¹ca œrodki znaczkina w kolenych klatkach dla video I
k=1;  
while ~isDone(I)
    next_frame1=step(I);
    for yyi=-zasieg:zasieg
        for xxi=-zasieg:zasieg
            com_block1=next_frame1(round(srodki1(k,2)-otoczenie1+yyi):round(srodki1(k,2)+otoczenie1+yyi),round(srodki1(k,1)-otoczenie1+xxi):round(srodki1(k,1)+otoczenie1+xxi),:);
              MSEtab(xxi+zasieg+1,yyi+zasieg+1) = sum(sum(sum((ref_block1-com_block1).^2)));
        end
    end
    tymczasowa = MSEtab(:);
    [minimum, index] = min(tymczasowa(:));
    [rzad, kolumna] = ind2sub(size(MSEtab),index);
  
    srodki1(k+1,1) = srodki1(k,1) + (rzad-1)- zasieg;
    srodki1(k+1,2) = srodki1(k,2) + (kolumna-1) - zasieg;

    clear MSEtab;
    k=k+1;    
end

%% 
%Pêtla wyznaczaj¹ca œrodki znaczkina w kolenych klatkach dla video II
l=1;  
while ~isDone(I2)
    next_frame2=step(I2);
    for yyi2=-zasieg:zasieg
        for xxi2=-zasieg:zasieg
            com_block2=next_frame2(round(srodki2(l,2)-otoczenie2+yyi2):round(srodki2(l,2)+otoczenie2+yyi2),round(srodki2(l,1)-otoczenie2+xxi2):round(srodki2(l,1)+otoczenie2+xxi2),:);
            MSEtab2(xxi2+zasieg+1,yyi2+zasieg+1) = sum(sum(sum((ref_block2-com_block2).^2)));
        end
    end
    tymczasowa2 = MSEtab2(:);
    [minimum2, index2] = min(tymczasowa2(:));
    [rzad2, kolumna2] = ind2sub(size(MSEtab2),index2);
  
    srodki2(l+1,1) = srodki2(l,1) + (rzad2-1)- zasieg;
    srodki2(l+1,2) = srodki2(l,2) + (kolumna2-1) - zasieg;
 
    clear MSEtab2;
    l=l+1;    
end

%%
figure
subplot(1,2,1) 
plot(srodki1(:,1),srodki1(:,2),'m');
set(gca, 'Ydir', 'reverse')
title('Trajektoria œrodka znacznika nagrania I');
xlabel('x');
ylabel('y');
axis([250 450 200 400])
subplot(1,2,2) 
plot(srodki2(:,1),srodki2(:,2));
set(gca, 'Ydir', 'reverse')
title('Trajektoria œrodka znacznika nagrania II');
xlabel('x');
ylabel('y');
axis([300 480 200 400])

%%
%nagrania z zaznaczonymi wspó³rzêdnymi markera (viedo I)
% https://uk.mathworks.com/help/vision/ref/inserttext.html?s_tid=srchtitle

release(I);
release(I2);
tekst = cell(nframes2,1);

for o=1:nframes2
    tekst{o} = [' X=' num2str(srodki2(l,1),'%0.2f') ' Y=' num2str(srodki2(l,2),'%0.2f')];
end
m=1;
while ~isDone(I2)
    frame=step(I2);
    out2 = insertText(frame,srodki2(m,:),tekst(m));
    step(videoPlayer,out2);
    m=m+1;
end

release(videoPlayer);
release(I);
release(I2);

%%
%nagranie ze znacznikiem (video II)
% https://uk.mathworks.com/help/vision/ref/insertmarker.html
release(videoPlayer);
release(I);
release(I2);

n=1;

while ~isDone(I)
    frame=step(I);
    out = insertMarker(frame,srodki1(n,:),'*', 'size',6);
    step(videoPlayer,out);
    n=n+1;
end

release(videoPlayer);
release(I);
release(I2);

%% 
%inrepolacja 
%xq - iloœæ punktów ktore chcemy z y2

if nframes1>nframes2
    xq = linspace(1, nframes2, nframes1);
    yy_Im2 = interp1(1:nframes2, srodki2(:,2), xq)';
    xx_Im2 = interp1(1:nframes2, srodki2(:,1), xq)';
elseif nframes1<nframes2
    xq = linspace(1, nframes1, nframes2);
    yy_Im1 = interp1(1:nframes1,srodki1(:,2),xq)';
    xx_Im1 = interp1(1:nframes1,srodki1(:,1),xq)';
end 

% tworzenie macierzy wspó³rzêdnych dla dwóch nagrañ

if nframes1>nframes2
    srodki2 = [xx_Im2, yy_Im2];
elseif nframes1<nframes2
    srodki1 = [xx_Im1, yy_Im1];
end 
%%
figure
plot3(srodki2(:,1), srodki2(:,2),srodki1(:,1))
set(gca, 'Ydir', 'reverse')
title('Trójwymiarowa trajektoria znacznika.')
xlabel('x');
ylabel('y');
zlabel('z');
axis([300 500 200 400 300 500]);
grid on

%%
%skalowanie 

Skalowanie_pion_bok = vision.VideoFileReader('skal_pion_bok.avi');
Skalowanie_poziom_bok = vision.VideoFileReader('skal_poziom_bok.avi');
Skalowanie_pion_front = vision.VideoFileReader('skal_pion_front.avi');
Skalowanie_poziom_front = vision.VideoFileReader('skal_poziom_front.avi');

S1 = Skalowanie_pion_bok();
S2 = Skalowanie_poziom_bok();
S3 = Skalowanie_pion_front();
S4 = Skalowanie_poziom_front();
 
figure
imshow(S1);
title('zmierz d³ugoœæ przedmiotu skaluj¹cego');
h = imdistline(gca);
api = iptgetapi(h);
pause();
dlugosc1 = api.getDistance();

imshow(S2);
title('zmierz d³ugoœæ przedmiotu skaluj¹cego');
h2 = imdistline(gca);
api = iptgetapi(h2);
pause();
dlugosc2 = api.getDistance();

imshow(S3);
title('zmierz d³ugoœæ przedmiotu skaluj¹cego');
h3 = imdistline(gca);
api = iptgetapi(h3);
pause();
dlugosc3 = api.getDistance();
 
imshow(S4);
title('zmierz d³ugoœæ przedmiotu skaluj¹cego');
h4 = imdistline(gca);
api = iptgetapi(h4);
pause();
dlugosc4 = api.getDistance();


%%
%skalowanie macierzy œrodki

srodki1_skal = srodki1.*znana_dlugosc./dlugosc3;
srodki2_skal = srodki2.*znana_dlugosc./dlugosc1;
     
figure
subplot(1,2,1) 
plot(srodki1_skal(:,1),srodki1_skal(:,2),'m');
set(gca, 'Ydir', 'reverse')
title('Trajektoria œrodka znacznika nagrania I');
xlabel('odleg³oœæ x [mm]');
ylabel('odlegloœæ y [mm]');
axis([800 1200 600 1000]);

subplot(1,2,2) 
plot(srodki2_skal(:,1),srodki2_skal(:,2));
set(gca, 'Ydir', 'reverse')
title('Trajektoria œrodka znacznika nagrania II');
xlabel('odleg³oœæ x [mm]');
ylabel('odlegloœæ y [mm]');
axis([850 1250 600 1000]);

figure
plot3(srodki2_skal(:,1), srodki2_skal(:,2),srodki1_skal(:,1))
set(gca, 'Ydir', 'reverse')
title('Trójwymiarowa trajektoria znacznika.')
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
axis([800 1200 600 1000 800 1200]);
grid on
%%
%dla nagrania z p³aszczyzny bocznej
%wzd³u¿ osi y
[min_y, index_miny] = min(srodki2(:,2));
[max_y, index_maxy] = max(srodki2(:,2));
%wzd³u¿ osi x
[min_x, index_minx] = min(srodki2(:,1));
[max_x, index_maxx] = max(srodki2(:,1));

figure
plot(srodki2(:,1), srodki2(:,2));
set(gca, 'Ydir', 'reverse')
title('Trajektoria œrodka znacznika nagrania I');
xlabel('x');
ylabel('y');

hold on
plot(srodki2(index_miny,1), min_y, 'o', srodki2(index_maxy,1), max_y, 'o')
hold on
plot(min_x, srodki2(index_minx,2), 'o', max_x, srodki2(index_maxx,2), 'o')

srednica_pion=(max_y-min_y)*znana_dlugosc/dlugosc1;
srednica_poziom=(max_x-min_x)*znana_dlugosc/dlugosc2;

fprintf('Odleg³oœæ pomiêdzy maksymaln¹ a minimaln¹ wartoœci¹ w osi y nagrania z p³aszczyzny bocznej wynosi %4.2f mm.\n', srednica_pion);
fprintf('Odleg³oœæ pomiêdzy maksymaln¹ a minimaln¹ wartoœci¹ w osi x nagrania z p³aszczyzny bocznej %4.2f mm.\n', srednica_poziom);

%%
%dla nagrania z p³aszczyzny frontowej
[min_y2, index_miny2] = min(srodki1(:,2));
[max_y2, index_maxy2] = max(srodki1(:,2));
%wzd³u¿ osi x
[min_x2, index_minx2] = min(srodki1(:,1));
[max_x2, index_maxx2] = max(srodki1(:,1));

figure
plot(srodki1(:,1), srodki1(:,2));
set(gca, 'Ydir', 'reverse')
title('Trajektoria œrodka znacznika nagrania II');
xlabel('x');
ylabel('y');

hold on
plot(srodki1(index_miny2,1), min_y2, 'o', srodki1(index_maxy2,1), max_y2, 'o')
hold on
plot(min_x2, srodki1(index_minx2,2), 'o', max_x2, srodki1(index_maxx2,2), 'o')

srednica_pion2=(max_y2-min_y2)*znana_dlugosc/dlugosc3;
srednica_poziom2=(max_x2-min_x2)*znana_dlugosc/dlugosc4;

fprintf('Odleg³oœæ pomiêdzy maksymaln¹ a minimaln¹ wartoœci¹ w osi y nagrania z p³aszczyzny frontowej wynosi %4.2f mm.\n', srednica_pion2);
fprintf('Odleg³oœæ pomiêdzy maksymaln¹ a minimaln¹ wartoœci¹ w osi x nagrania z p³aszczyzny frontowej %4.2f mm.\n', srednica_poziom2);