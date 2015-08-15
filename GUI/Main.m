function varargout = Main(varargin)
% MAIN MATLAB code for Main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main

% Last Modified by GUIDE v2.5 17-Feb-2015 17:25:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_OutputFcn, ...
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


% --- Executes just before Main is made visible.
function Main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main (see VARARGIN)

% Choose default command line output for Main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes Main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;

global setup;
global soundfield;
global quiet;
global bright;

Res          = str2num(get(handles.numRes, 'String'));
pw_angle     = get(handles.numPWangle,'Value');
f            = str2num(get(handles.numFrequency,'String'));
QuietRadius  = str2num(get(handles.numQuietRadius,'String'));
BrightRadius = str2num(get(handles.numBrightRadius,'String'));

quiet  = Orthogonal_Basis_Expansion.spatial_zone(f, 0, QuietRadius, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(f, 0, BrightRadius, 'pw', 1.0, pw_angle);
quiet.res  = Res;
bright.res = quiet.res;
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');
bright = bright.setDesiredSoundfield(true);

%%
BrightLocAngle = str2num(get(handles.numBrightLocAngle,'String'));
BrightLocDist  = str2num(get(handles.numBrightLocDist,'String'));
QuietLocAngle  = str2num(get(handles.numQuietLocAngle,'String'));
QuietLocDist   = str2num(get(handles.numQuietLocDist,'String'));

soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  QuietLocDist, QuietLocAngle);
soundfield = soundfield.addSpatialZone(bright, BrightLocDist, BrightLocAngle);

%%
reprodRad    = str2num(get(handles.numRadius,'String'));
N            = str2num(get(handles.numN,'String'));
BrightWeight = str2num(get(handles.numBrightWeight,'String'));
QuietWeight  = str2num(get(handles.numQuietWeight,'String'));
OtherWeight  = str2num(get(handles.numOtherWeight,'String'));

soundfield.BrightZ_Weight     = BrightWeight;
soundfield.QuietZ_Weight      = QuietWeight;
soundfield.UnattendedZ_Weight = OtherWeight;

soundfield = soundfield.setN(N);
soundfield = soundfield.createSoundfield('DEBUG', reprodRad);

%%
N_loudspkrs    = str2num(get(handles.numLoudspeakers,'String'));
Rad_loudspkrs  = str2num(get(handles.numLoudspeakerRadius,'String'));
Arc_loudspkrs  = str2num(get(handles.numLoudspeakerArc,'String'));
First_loudspkr = str2num(get(handles.numLoudspeakerFirst,'String'));

setup = Speaker_Setup.loudspeaker_setup;
setup = setup.addMultizone_Soundfield(soundfield);
setup.Loudspeaker_Count = N_loudspkrs;
setup.Speaker_Arc_Angle = Arc_loudspkrs;
setup.Angle_FirstSpeaker = First_loudspkr;
setup = setup.setRadius(Rad_loudspkrs);

setup = setup.calc_Loudspeaker_Weights();
setup = setup.reproduceSoundfield('DEBUG');

%%
a = axes('Tag', 'axes1');
axes(handles.axes1);
delete(a)
global h;
h = setup.plotSoundfield();




% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on slider movement. 
function numPWangle_Callback(hObject, eventdata, handles)
% hObject    handle to numPWangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.txtPWAngle,'String',[num2str(get(handles.numPWangle,'Value')) '°']);
pushbutton1_Callback(handles.pushbutton1, eventdata, handles);


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function numPWangle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numPWangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global animate_state;
animate_state = false;
set(handles.pushbutton2, 'Enable', 'off');
set(handles.pushbutton3, 'Enable', 'on');
set(handles.pushbutton1, 'Enable', 'on');
set(handles.numPWangle, 'Enable', 'on');
set(handles.numCalcMinReqSpkrs, 'Enable', 'on');



function numRes_Callback(hObject, eventdata, handles)
% hObject    handle to numRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numRes as text
%        str2double(get(hObject,'String')) returns contents of numRes as a double


% --- Executes during object creation, after setting all properties.
function numRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)phase = 0;
global setup;
global h;
global animate_state;
animate_state=true;

set(handles.pushbutton2, 'Enable', 'on');
set(handles.pushbutton3, 'Enable', 'off');
set(handles.pushbutton1, 'Enable', 'off');
set(handles.numPWangle, 'Enable', 'off');
set(handles.numCalcMinReqSpkrs, 'Enable', 'off');

frames = [];
for phase = 1:1:36
    frames(:,:,phase)=setup.Soundfield_reproduced * exp(-1i*phase/18*pi);
end
phase = 0;
while (animate_state)
    phase = phase + 1;
    if phase == 37
        phase = 1;
    end
    set(h,'ZData', frames(:,:,phase) );
    drawnow();
end



function numFrequency_Callback(hObject, eventdata, handles)
% hObject    handle to numFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numFrequency as text
%        str2double(get(hObject,'String')) returns contents of numFrequency as a double


% --- Executes during object creation, after setting all properties.
function numFrequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numRadius_Callback(hObject, eventdata, handles)
% hObject    handle to numRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numRadius as text
%        str2double(get(hObject,'String')) returns contents of numRadius as a double


% --- Executes during object creation, after setting all properties.
function numRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numN_Callback(hObject, eventdata, handles)
% hObject    handle to numN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numN as text
%        str2double(get(hObject,'String')) returns contents of numN as a double


% --- Executes during object creation, after setting all properties.
function numN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numBrightWeight_Callback(hObject, eventdata, handles)
% hObject    handle to numBrightWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numBrightWeight as text
%        str2double(get(hObject,'String')) returns contents of numBrightWeight as a double


% --- Executes during object creation, after setting all properties.
function numBrightWeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numBrightWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numQuietWeight_Callback(hObject, eventdata, handles)
% hObject    handle to numQuietWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numQuietWeight as text
%        str2double(get(hObject,'String')) returns contents of numQuietWeight as a double


% --- Executes during object creation, after setting all properties.
function numQuietWeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numQuietWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numOtherWeight_Callback(hObject, eventdata, handles)
% hObject    handle to numOtherWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numOtherWeight as text
%        str2double(get(hObject,'String')) returns contents of numOtherWeight as a double


% --- Executes during object creation, after setting all properties.
function numOtherWeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numOtherWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numBrightRadius_Callback(hObject, eventdata, handles)
% hObject    handle to numBrightRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numBrightRadius as text
%        str2double(get(hObject,'String')) returns contents of numBrightRadius as a double


% --- Executes during object creation, after setting all properties.
function numBrightRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numBrightRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numBrightLocAngle_Callback(hObject, eventdata, handles)
% hObject    handle to numBrightLocAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numBrightLocAngle as text
%        str2double(get(hObject,'String')) returns contents of numBrightLocAngle as a double


% --- Executes during object creation, after setting all properties.
function numBrightLocAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numBrightLocAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numBrightLocDist_Callback(hObject, eventdata, handles)
% hObject    handle to numBrightLocDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numBrightLocDist as text
%        str2double(get(hObject,'String')) returns contents of numBrightLocDist as a double


% --- Executes during object creation, after setting all properties.
function numBrightLocDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numBrightLocDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numQuietRadius_Callback(hObject, eventdata, handles)
% hObject    handle to numQuietRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numQuietRadius as text
%        str2double(get(hObject,'String')) returns contents of numQuietRadius as a double


% --- Executes during object creation, after setting all properties.
function numQuietRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numQuietRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numQuietLocAngle_Callback(hObject, eventdata, handles)
% hObject    handle to numQuietLocAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numQuietLocAngle as text
%        str2double(get(hObject,'String')) returns contents of numQuietLocAngle as a double


% --- Executes during object creation, after setting all properties.
function numQuietLocAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numQuietLocAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numQuietLocDist_Callback(hObject, eventdata, handles)
% hObject    handle to numQuietLocDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numQuietLocDist as text
%        str2double(get(hObject,'String')) returns contents of numQuietLocDist as a double


% --- Executes during object creation, after setting all properties.
function numQuietLocDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numQuietLocDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numLoudspeakers_Callback(hObject, eventdata, handles)
% hObject    handle to numLoudspeakers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numLoudspeakers as text
%        str2double(get(hObject,'String')) returns contents of numLoudspeakers as a double


% --- Executes during object creation, after setting all properties.
function numLoudspeakers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numLoudspeakers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numLoudspeakerRadius_Callback(hObject, eventdata, handles)
% hObject    handle to numLoudspeakerRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numLoudspeakerRadius as text
%        str2double(get(hObject,'String')) returns contents of numLoudspeakerRadius as a double


% --- Executes during object creation, after setting all properties.
function numLoudspeakerRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numLoudspeakerRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numLoudspeakerArc_Callback(hObject, eventdata, handles)
% hObject    handle to numLoudspeakerArc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numLoudspeakerArc as text
%        str2double(get(hObject,'String')) returns contents of numLoudspeakerArc as a double


% --- Executes during object creation, after setting all properties.
function numLoudspeakerArc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numLoudspeakerArc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numLoudspeakerFirst_Callback(hObject, eventdata, handles)
% hObject    handle to numLoudspeakerFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numLoudspeakerFirst as text
%        str2double(get(hObject,'String')) returns contents of numLoudspeakerFirst as a double


% --- Executes during object creation, after setting all properties.
function numLoudspeakerFirst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numLoudspeakerFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in numCalcMinReqSpkrs.
function numCalcMinReqSpkrs_Callback(hObject, eventdata, handles)
% hObject    handle to numCalcMinReqSpkrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Radius   = str2num(get(handles.numRadius,'String'));
f_max    = str2num(get(handles.numMaxFrequency,'String'));
arc      = str2num(get(handles.numLoudspeakerArc,'String'));

set(handles.numLoudspeakers , 'String', num2str(ceil((2*ceil(f_max/343 * exp(1) * Radius / 2) + 1) * arc/360)));



function numMaxFrequency_Callback(hObject, eventdata, handles)
% hObject    handle to numMaxFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numMaxFrequency as text
%        str2double(get(hObject,'String')) returns contents of numMaxFrequency as a double


% --- Executes during object creation, after setting all properties.
function numMaxFrequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numMaxFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
