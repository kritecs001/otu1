// Execute from Scilab. Results are written next to this script.
function data=simulateLab(diagram)
    info=scicos_simulate(diagram,list(),struct(),"nw");
    if ~exists("lab_data","local") then
        error("Xcos did not return lab_data.");
    end
    data=[lab_data.time, lab_data.values];
    if size(data,1)<2 | or(isnan(data)) | or(isinf(data)) then
        error("Invalid Xcos output.");
    end
endfunction

function modeDiagram=setLabMode(diagram,mode)
    modeDiagram=diagram;
    context=diagram.props.context;
    context(1)="mode="+string(mode);
    modeDiagram.props.context=context;
endfunction

function writeTimeData(root,out)
    names=["01_integrator_k";"02_aperiodic_k";"03_aperiodic_T"; ..
           "04_oscillatory_k";"05_oscillatory_T";"06_oscillatory_d"; ..
           "07_differentiator";"09_derivative_comparison"];
    for j=1:size(names,1)
        name=names(j);
        diagram=xcosDiagramToScilab(root+name+".xcos");
        for mode=0:1
            data=simulateLab(setLabMode(diagram,mode));
            if mode==0 then suffix="_step"; else suffix="_pulse"; end
            csvWrite(data,out+name+suffix+".csv",",",".","%.12g");
        end
        mprintf("%s OK\n",name);
    end
endfunction

function writeFreqData(root,out)
    model=xcosDiagramToScilab(root+"08_sine_frequency.xcos");
    frequencies=[0.1;0.2;0.5;1;2;10/3;5;10;20;50];
    freq=zeros(size(frequencies,1),9);
    for j=1:length(frequencies)
        w=frequencies(j);
        P=2*%pi/w;
        dt=P/160;
        n=max(8,ceil(12*0.3/P)+5);
        model.props.tf=n*P+dt/10;
        model.props.context=[msprintf("omega=%.17g",w);msprintf("dt=%.17g",dt)];
        data=simulateLab(model);
        start=find(data(:,1)>=(n-4)*P-dt/10);
        t=data(start,1);
        // The offset term removes the startup constant of the ideal integrator.
        B=[ones(t),sin(w*t),cos(w*t)];
        c=B\data(start,2:4);
        H=(c(2,2:3)+%i*c(3,2:3))/(c(2,1)+%i*c(3,1));
        exact=[5/(%i*w),5*%i*w/(1+0.3*%i*w)];
        freq(j,:)=[w,abs(H(1)),atan(imag(H(1)),real(H(1)))*180/%pi, ..
                   abs(exact(1)),atan(imag(exact(1)),real(exact(1)))*180/%pi, ..
                   20*log10(abs(H(2))),atan(imag(H(2)),real(H(2)))*180/%pi, ..
                   20*log10(abs(exact(2))),atan(imag(exact(2)),real(exact(2)))*180/%pi];
    end
    csvWrite(freq,out+"frequency.csv",",",".","%.12g");
    mprintf("08_sine_frequency OK\n");
endfunction

root=get_absolute_file_path("run_all.sce");
out=root+"results";
if ~isdir(out) && ~createdir(out) then
    error("Could not create results directory: "+out);
end
out=out+"/";
// Initialize the optional Modelica scratch directory without a shell call.
modelicaCache=TMPDIR+"/modelica";
if ~isdir(modelicaCache) then createdir(modelicaCache); end
loadXcosLibs();
writeTimeData(root,out);
writeFreqData(root,out);
mprintf("CSV: %s\n",out);
