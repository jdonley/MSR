function res = helloworld
%res = 'Hello World!';
h=figure('Visible','off');
text(0.2,0.5,'Hello World!','fontsize',40);
res = webfigure(h);
end

