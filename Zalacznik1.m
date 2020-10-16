%PROGRAM DO WYZNACZANIA TRÓJWYMIAROWEJ TRAJEKTORII ZNACZNIKA NA CIELE
%Aleksandra Michalska

function varargout = Zalacznik1(varargin)
% ZALACZNIK1 MATLAB code for Zalacznik1.fig
%      ZALACZNIK1, by itself, creates a new ZALACZNIK1 or raises the existing
%      singleton*.
%
%      H = ZALACZNIK1 returns the handle to a new ZALACZNIK1 or the handle to
%      the existing singleton*.
%
%      ZALACZNIK1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZALACZNIK1.M with the given input arguments.
%
%      ZALACZNIK1('Property','Value',...) creates a new ZALACZNIK1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Zalacznik1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Zalacznik1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Zalacznik1

% Last Modified by GUIDE v2.5 26-Apr-2019 13:41:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Zalacznik1_OpeningFcn, ...
                   'gui_OutputFcn',  @Zalacznik1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Zalacznik1 is made visible.
function Zalacznik1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Zalacznik1 (see VARARGIN)

% Choose default command line output for Zalacznik1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.zatwierdz_nazwy_skal,'enable','off')
set(handles.zatwierdz_skal,'enable','off')
set(handles.nagranie_bok,'enable','off')
set(handles.nagranie_front,'enable','off')
set(handles.wartosc_zasieg,'enable','off')
set(handles.zatwierdz_nazwy_nagran,'enable','off')
set(handles.zatwierdz_pierwsze_klatki,'enable','off')
set(handles.mat_plik,'enable','off')

% UIWAIT makes Zalacznik1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Zalacznik1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in wybierz_pion_bok.
function wybierz_pion_bok_Callback(hObject, eventdata, handles)
% hObject    handle to wybierz_pion_bok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nskal_pion_bok, path1] = uigetfile('*.avi');           %wczytuje nagranie
set(handles.nazwa_pion_bok, 'string', nskal_pion_bok);  %wyswietla nazwe nagrania
assignin('base', 'skal_pion_bok', nskal_pion_bok);      %zapisuje do workspace

% --- Executes on button press in wybierz_poziom_bok.
function wybierz_poziom_bok_Callback(hObject, eventdata, handles)
% hObject    handle to wybierz_poziom_bok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nskal_poziom_bok, path2] = uigetfile('*.avi');
set(handles.nazwa_poziom_bok, 'string', nskal_poziom_bok);
assignin('base', 'skal_poziom_bok', nskal_poziom_bok);

% --- Executes on button press in wybierz_pion_front.
function wybierz_pion_front_Callback(hObject, eventdata, handles)
% hObject    handle to wybierz_pion_front (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nskal_pion_front, path3] = uigetfile('*.avi');
set(handles.nazwa_pion_front, 'string', nskal_pion_front);
assignin('base', 'skal_pion_front', nskal_pion_front);

% --- Executes on button press in wybierz_poziom_front.
function wybierz_poziom_front_Callback(hObject, eventdata, handles)
% hObject    handle to wybierz_poziom_front (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nskal_poziom_front, path4] = uigetfile('*.avi');
set(handles.nazwa_poziom_front, 'string', nskal_poziom_front);
assignin('base', 'skal_poziom_front', nskal_poziom_front);

function dlugosc_Callback(hObject, eventdata, handles)
% hObject    handle to dlugosc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dlugosc as text
%        str2double(get(hObject,'String')) returns contents of dlugosc as a double
znana_dlugosc = str2double(get(handles.dlugosc, 'string'));
set(handles.dlugosc2, 'string', znana_dlugosc);
save znana_dlugosc znana_dlugosc;

set(handles.zatwierdz_nazwy_skal,'enable','on')

% --- Executes during object creation, after setting all properties.
function dlugosc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dlugosc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in zatwierdz_nazwy_skal.
function zatwierdz_nazwy_skal_Callback(hObject, eventdata, handles)
% hObject    handle to zatwierdz_nazwy_skal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

B1 = get(handles.nazwa_pion_bok, 'String');
B2 = get(handles.nazwa_poziom_bok, 'String');
B3 = get(handles.nazwa_pion_front, 'String');
B4 = get(handles.nazwa_poziom_front, 'String');

if isempty(B1)
    errordlg('Wgraj nagranie z p³aszczyzny bocznej, które skaluje w pionie.','B³¹d')
    return
elseif isempty(B2)
    errordlg('Wgraj nagranie z p³aszczyzny bocznej, które skaluje w poziomie.','B³¹d')
    return
elseif isempty(B3)
    errordlg('Wgraj nagranie z p³aszczyzny frontowej, które skaluje w pionie.','B³¹d')
    return
elseif isempty(B4)
    errordlg('Wgraj nagranie z p³aszczyzny frontowej, które skaluje w poziomie.','B³¹d')
    return
else
    set(handles.text_skal, 'string', 'Wska¿ d³ugoœæ przedmiotu skaluj¹cego na poni¿szych obrazach.');

    nskal_pion_bok = get(handles.nazwa_pion_bok, 'string');
    Skalowanie_pion_bok = vision.VideoFileReader(nskal_pion_bok);
    S1 = Skalowanie_pion_bok();
    axes(handles.skal_pion_bok)
    imshow(S1);
    h = imdistline(gca);
    api = iptgetapi(h);
    pause();
    dlugosc1 = api.getDistance();

    nskal_poziom_bok = get(handles.nazwa_poziom_bok, 'string');
    Skalowanie_poziom_bok = vision.VideoFileReader(nskal_poziom_bok);
    S2 = Skalowanie_poziom_bok();
    axes(handles.skal_poziom_bok)
    imshow(S2);
    h2 = imdistline(gca);
    api = iptgetapi(h2);
    pause();
    dlugosc2 = api.getDistance();

    nskal_pion_front = get(handles.nazwa_pion_front, 'string');
    Skalowanie_pion_front = vision.VideoFileReader(nskal_pion_front);
    S3 = Skalowanie_pion_front();
    axes(handles.skal_pion_front)
    imshow(S3);
    h3 = imdistline(gca);
    api = iptgetapi(h3);
    pause();
    dlugosc3 = api.getDistance();

    nskal_poziom_front = get(handles.nazwa_poziom_front, 'string');
    Skalowanie_poziom_front = vision.VideoFileReader(nskal_poziom_front);
    S4 = Skalowanie_poziom_front();
    axes(handles.skal_poziom_front)
    imshow(S4);
    h4 = imdistline(gca);
    api = iptgetapi(h4);
    pause();
    dlugosc4 = api.getDistance();

    save dlugosci dlugosc1 dlugosc2 dlugosc3 dlugosc4;
    set(handles.zatwierdz_skal,'enable','on')
end

% --- Executes on button press in zatwierdz_skal.
function zatwierdz_skal_Callback(hObject, eventdata, handles)
% hObject    handle to zatwierdz_skal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.nagranie_bok,'enable','on')
set(handles.nagranie_front,'enable','on')
set(handles.wartosc_zasieg,'enable','on')


% --- Executes on button press in robocza.
function robocza_Callback(hObject, eventdata, handles)
% hObject    handle to robocza (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in mat_plik.
function mat_plik_Callback(hObject, eventdata, handles)
% hObject    handle to mat_plik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

wykres = get(handles.wykres3D);
save wykres3D wykres;
set(handles.zapis_mat, 'String', 'Zapisano');

% --- Executes on button press in zatwierdz_pierwsze_klatki.
function zatwierdz_pierwsze_klatki_Callback(hObject, eventdata, handles)
% hObject    handle to zatwierdz_pierwsze_klatki (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load zmienne1;
load dlugosci;
load zasieg;
load znana_dlugosc;

srodki1=zeros(nframes1,2);
srodki2=zeros(nframes2,2);

srodki1(1,1)=base_xp1;
srodki1(1,2)=base_yp1;
srodki2(1,1)=base_xp2;
srodki2(1,2)=base_yp2;

nagranieBok = get(handles.nazwa_bok, 'string');
I = vision.VideoFileReader(nagranieBok);
base_frame1 = I();
ref_block1=base_frame1(round(srodki1(1,2)-kernel1):round(srodki1(1,2)+kernel1),round(srodki1(1,1)-kernel1):round(srodki1(1,1)+kernel1),:);

nagranieFront = get(handles.nazwa_front, 'string');
I2 = vision.VideoFileReader(nagranieFront);
base_frame2 = I2(); 
ref_block2=base_frame2(round(srodki2(1,2)-kernel2):round(srodki2(1,2)+kernel2),round(srodki2(1,1)-kernel2):round(srodki2(1,1)+kernel2),:);


%Pêtla wyznaczaj¹ca œrodki znaczkina w kolenych klatkach dla video I
k=1;  
while ~isDone(I)
    next_frame1=step(I);
    for yyi=-zasieg:zasieg
        for xxi=-zasieg:zasieg
            com_block1=next_frame1(round(srodki1(k,2)-kernel1+yyi):round(srodki1(k,2)+kernel1+yyi),round(srodki1(k,1)-kernel1+xxi):round(srodki1(k,1)+kernel1+xxi),:);
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

%Pêtla wyznaczaj¹ca œrodki znaczkina w kolenych klatkach dla video II
l=1;  
while ~isDone(I2)
    next_frame2=step(I2);
    for yyi2=-zasieg:zasieg
        for xxi2=-zasieg:zasieg
            com_block2=next_frame2(round(srodki2(l,2)-kernel2+yyi2):round(srodki2(l,2)+kernel2+yyi2),round(srodki2(l,1)-kernel2+xxi2):round(srodki2(l,1)+kernel2+xxi2),:);
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

if nframes1>nframes2
    srodki2 = [xx_Im2, yy_Im2];
elseif nframes1<nframes2
    srodki1 = [xx_Im1, yy_Im1];
end 

srodki1_skal = srodki1.*znana_dlugosc./dlugosc3;

srodki2_skal = srodki2.*znana_dlugosc./dlugosc1;

set(handles.mat_plik,'enable','on')


axes(handles.wykres3D);
plot3(srodki2_skal(:,1), srodki2_skal(:,2),srodki1_skal(:,1))
set(gca, 'Ydir', 'reverse')
title('Trójwymiarowa trajektoria znacznika.')
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
axis([800 1200 450 850 800 1200]);


% --- Executes on button press in nagranie_bok.
function nagranie_bok_Callback(hObject, eventdata, handles)
% hObject    handle to nagranie_bok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nagranieBok, path5] = uigetfile('*.avi');
set(handles.nazwa_bok, 'string', nagranieBok);
assignin('base', 'nagranie_bok', nagranieBok);

% --- Executes on button press in nagranie_front.
function nagranie_front_Callback(hObject, eventdata, handles)
% hObject    handle to nagranie_front (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nagranieFront, path6] = uigetfile('*.avi');
set(handles.nazwa_front, 'string', nagranieFront);
assignin('base', 'nagranie_front', nagranieFront);

function wartosc_zasieg_Callback(hObject, eventdata, handles)
% hObject    handle to wartosc_zasieg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wartosc_zasieg as text
%        str2double(get(hObject,'String')) returns contents of wartosc_zasieg as a double
zasieg = str2double(get(handles.wartosc_zasieg, 'string'));
set(handles.zasieg2, 'string', zasieg);
save zasieg zasieg;

set(handles.zatwierdz_nazwy_nagran,'enable','on')



% --- Executes during object creation, after setting all properties.
function wartosc_zasieg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wartosc_zasieg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in zatwierdz_nazwy_nagran.
function zatwierdz_nazwy_nagran_Callback(hObject, eventdata, handles)
% hObject    handle to zatwierdz_nazwy_nagran (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


B5 = get(handles.nazwa_bok, 'String');
B6 = get(handles.nazwa_front, 'String');

if isempty(B5)
    errordlg('Wgraj nagranie z p³aszczyzny bocznej.','B³¹d')
    return
elseif isempty(B6)
    errordlg('Wgraj nagranie z p³aszczyzny frontowej.','B³¹d')
    return
else
    nagranieBok = get(handles.nazwa_bok, 'string');
    I = vision.VideoFileReader(nagranieBok);
    Im1 = VideoReader(nagranieBok);
    nframes1 = Im1.NumberOfFrames;

    base_frame1 = I();   
    axes(handles.pierwsza_bok);
    imshow(base_frame1);
    title({'Pierwsza klatka nagrania I.';'1. Zaznacz obszar wokó³ znacznika.';'2. Zaznacz œrodek markera.'});

    objectRegion=round(getPosition(imrect));
    kernel1=objectRegion(1,3)/2;
    [base_xp1,base_yp1]=ginput(1);

    nagranieFront = get(handles.nazwa_front, 'string');
    I2 = vision.VideoFileReader(nagranieFront);
    Im2 = VideoReader(nagranieFront);
    nframes2 = Im2.NumberOfFrames;

    base_frame2 = I2();   
    axes(handles.pierwsza_front);
    imshow(base_frame2);
    title({'Pierwsza klatka nagrania II.';'1. Zaznacz obszar wokó³ znacznika.';'2. Zaznacz œrodek markera.'});

    objectRegion2=round(getPosition(imrect));
    kernel2=objectRegion2(1,3)/2;
    [base_xp2,base_yp2]=ginput(1);

    save zmienne1 kernel1 base_xp1 nframes1 base_yp1 kernel2 base_xp2 base_yp2 nframes2;

    set(handles.zatwierdz_pierwsze_klatki,'enable','on')
end
