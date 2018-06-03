frames = 100;  %Zoveel frames
y = sin(linspace(0,10,100)); %Functie die geanimeerd wordt
h = figure;
%set(h, 'Visible', 'off'); %Plot onzichtbaar, lijkt 'getframe' wel
%langzamer te maken
plot(y); %Plot begintoestand/initialiseer de assen
axis tight manual %Houdt dezelfde assen in het geval van 'hold on'
ax = gca; 
ax.NextPlot = 'replaceChildren'; %Houdt dezelfde assen nu ook bij vervolgplots

F(frames) = struct('cdata',[],'colormap',[]); %Structuur waarin de plot wordt opgeslagen


for j = 1:frames 
    X = cos(2*pi*j/frames)*y; %Veranderende functie 
    plot(X,'ro'); %Plot
    %drawnow %Als je het nu direct geplot wil
    
    F(j) = getframe(gcf); %Save frame
end


%fig = figure;
%movie(fig,F,2)

v = VideoWriter('testVideo.avi'); %Maak een video-file
open(v)
writeVideo(v,F) %Sla de frames op in de video
close(v)