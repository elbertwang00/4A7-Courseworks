function airload(filename)
% process the values in a file produced by airfoil.
% Assumes Xupper, CPupper, Xlower, CPlower, em, cl and cdv
% exist as globals. Values are returned into these variables.
% It assumes that 
% "Upper surface CP values" starts a table of values which is
% ended by "Lower surface CP values", which also starts another
% table ended by a blank line.
% Also assumes that em, cl and cdv are set at the start of a line.
%
% Sample usage: 
%   global Xupper CPupper Xlower CPlower em cl cdv
%   holgar Autorun.BRF.invis

global Xupper CPupper Xlower CPlower em cl cm cdv cd1 cd2 HTE Xsh Msh
if nargin < 1
  warning('You need to provide a filename');
  return
end


fid=fopen(filename);
if fid == -1
  warning('File does not exist');
  return
end

while 1
   line = fgetl(fid);
   if ~ischar(line)
     break
   end
   if strncmp('EM =',dblnk(line),4)
      strs=sscanf(line,'%s', 3);
      em=str2double(strs(4:end));
   end 
   if strncmp('CL =',dblnk(line),4)
      strs=sscanf(line,'%s', 6);
      cl=str2double(strs(4:10));
      cm=str2double(strs(14:end));
   end
   if strncmp('CDV =',dblnk(line),5)
      strs=sscanf(line,'%s', 3);
      cdv=str2double(strs(5:end));
   end
   if strncmp('CDV+CD1 =',dblnk(line),9)
      strs=sscanf(line,'%s', 6);
      cd1 = str2double(strs(9:15)) - cdv;
      cd2 = str2double(strs(24:end)) - cdv;
   end
   
   if strcmp('Upper surface CP values',dblnk(line))
      line = fgetl(fid);
      Xupper=[];
      CPupper=[];
      while 1
         line=fgetl(fid);
         if strcmp('Lower surface CP values',dblnk(line))
            break
         end
         strs=sscanf(line,'%f %f');
         Xupper=[Xupper; strs(1)];
         CPupper=[CPupper; strs(2)];
      end
      Xlower=[];
      CPlower=[];
      line=fgetl(fid);
      while 1
         line=fgetl(fid);
         if ~isstr(line)
	   break
	 end
         if strcmp('',line)
	   break
	 
	 end
	 if length(line) <2
	     break
	 end
         [strs,count]=sscanf(line,'%f %f');
         if (count == 2)
            Xlower=[Xlower; strs(1)];
            CPlower=[CPlower; strs(2)];
         end
      end 

   end 
end
fclose(fid);

FULlines = readlines('Autorun.FUL');

if strcmp(FULlines(end-1), " ITERATIVE PROCEDURE HAS DIVERGED")
    [em, cl, cm, cdv, cd1, cd2, HTE, Xsh, Msh] = deal(-1,-1,-1,-1,-0.0001,-0.0001,-1, -1, -1);
else
    HbarSection = find(strcmp(FULlines,'                                   UPPER-SURFACE BOUNDARY-LAYER DATA'));
    
    HTE = -1;

    for ii=HbarSection:1:length(FULlines)
        TEline = split(FULlines(ii));
        
        if length(TEline) == 7
            if str2double(TEline(2)) == 1
                HTE = str2double(TEline(3));
                break
            end
        end
    end

    XshLine = find(strcmp(FULlines, ' LOWER SURFACE.'));
    if ~isempty(XshLine)
        shockLine = split(FULlines(XshLine-4));
    
        if ~strcmp(shockLine(2), 'X(SHOCK)=')
            Xsh = -1;
            Msh = -1;
        else
            Xsh = str2double(shockLine(3));
            cpsh = CPupper(Xupper==Xsh);
            Msh = sqrt(2/0.4*((1+0.2*em^2)/(1+1.4*cpsh*em^2/2)^(0.4/1.4) - 1));
        end
    end
end



