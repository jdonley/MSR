function d = Dependencies
%DEPENDENCIES A structure array of all dependencies and associated info
d.Name=[];
d.URL=[];
d.Required=[];
d.IncludeSubDirs=[];

d(end+1).Name           = 'amtoolbox';
  d(end).URL            = 'http://amtoolbox.sourceforge.net/';
  d(end).Required       = true;
  d(end).IncludeSubDirs = true;
              
d(end+1).Name           = 'arrow';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/278-arrow';
              
d(end+1).Name           = 'cprintf';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window';
  d(end).Required       = true;
  
d(end+1).Name           = 'DSP_Tools';
  d(end).URL            = '';
  d(end).IncludeSubDirs = true;
  
d(end+1).Name           = 'export_fig';
  d(end).URL            = 'https://github.com/altmany/export_fig';
  d(end).Required       = true;
              
d(end+1).Name           = 'fastISM';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/25965-fast-simulation-of-acoustic-room-impulse-responses--image-source-method-';
              
d(end+1).Name           = 'GenerateFunctionMFile';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/27225-generate-a-new-function-m-file--with-documentation';
              
d(end+1).Name           = 'gridLegend';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/29248-gridlegend-a-multi-column-format-for-legends';
  d(end).Required       = true;
              
d(end+1).Name           = 'ISO226';
  d(end).URL            = 'http://au.mathworks.com/matlabcentral/fileexchange/7028-iso-226-equal-loudness-level-contour-signal';
              
d(end+1).Name           = 'm2html';
  d(end).URL            = 'https://www.artefact.tk/software/matlab/m2html/';
  d(end).Required       = true;
              
d(end+1).Name           = 'MCRoomSim';
  d(end).URL            = 'http://www.ee.usyd.edu.au/carlab/mcroomsim.htm';
              
d(end+1).Name           = 'mmstream2';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/38860-improved-2-d-streamlines';
              
d(end+1).Name           = 'mmx';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/37515-mmx-multithreaded-matrix-operations-on-n-d-matrices';
              
d(end+1).Name           = 'mtimesx';
  d(end).URL            = 'http://au.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support';
              
d(end+1).Name           = 'newFunction';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/27132-automatic-template-for-new-functions/';
              
d(end+1).Name           = 'Parfor_Progress';
  d(end).URL            = 'http://au.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor/';
  d(end).Required       = true;
              
d(end+1).Name           = 'PhaseUnwrapping2D';
  d(end).URL            = 'http://au.mathworks.com/matlabcentral/fileexchange/29497-goldsteinunwrap2d-r1';
              
d(end+1).Name           = 'playrec-master';
  d(end).URL            = 'https://github.com/PlayrecForMatlab/playrec';
  d(end).Required       = true;
              
d(end+1).Name           = 'python';
  d(end).URL            = '';
              
d(end+1).Name           = 'RIR_Generator';
  d(end).URL            = 'https://github.com/ehabets/RIR-Generator';
  d(end).Required       = true;
              
d(end+1).Name           = 'roomsim';
  d(end).URL            = 'https://sourceforge.net/projects/roomsim/';
              
d(end+1).Name           = 'scientific_colormaps';
  d(end).URL            = 'http://au.mathworks.com/matlabcentral/fileexchange/42450-custom-colormap';
  d(end).Required       = true;
              
d(end+1).Name           = 'Speech_Recognition_API';
  d(end).URL            = 'https://github.com/Uberi/speech_recognition';
              
d(end+1).Name           = 'SpeechSP_Tools';
  d(end).URL            = 'https://github.com/JacobD10/SoundZone_Tools';
  d(end).IncludeSubDirs = true;
              
d(end+1).Name           = 'sqdistance';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/32994-the-mincentropy-algorithm-for-alternative-clustering/content/minCEntropy/sqdistance.m';
              
d(end+1).Name           = 'SweptSineAnalysis';
  d(end).URL            = 'http://au.mathworks.com/matlabcentral/fileexchange/29187-swept-sine-analysis';
              
d(end+1).Name           = 'tightfigadv';
  d(end).URL            = 'https://au.mathworks.com/matlabcentral/fileexchange/65004-jacobd10-tightfigadv';
  d(end).Required       = true;
              
d(end+1).Name           = 'tightPlots';
  d(end).URL            = 'http://au.mathworks.com/matlabcentral/fileexchange/45380-tightplots';
  d(end).Required       = true;
              
d(end+1).Name           = 'voicebox';
  d(end).URL            = 'http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html';
  d(end).Required       = true;
    
  
d(1)=[];  
end

