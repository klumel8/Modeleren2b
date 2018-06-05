frames = 100;  %Zoveel frames
y = sin(linspace(0,10,100)); %Functie die geanimeerd wordt
y2 = sin(linspace(0,10,100))+2*sin(0.5*linspace(0,10,100));
subplot(1,2,1)
%plot(y); %Plot begintoestand/initialiseer de assen
axis tight manual %Houdt dezelfde assen in het geval van 'hold on'
ax = gca; 
ax.NextPlot = 'replaceChildren'; %Houdt dezelfde assen nu ook bij vervolgplots


subplot(1,2,2)
plot(y2)
axis tight manual %Houdt dezelfde assen in het geval van 'hold on'
ax = gca; 
ax.NextPlot = 'replaceChildren'; %Houdt dezelfde assen nu ook bij vervolgplots

F(frames) = struct('cdata',[],'colormap',[]); %Structuur waarin de plot wordt opgeslagen


for j = 1:frames 
    X = cos(2*pi*j/frames)*y; %Veranderende functie 
    X2 = cos(2*pi*j/frames)*y2;
    subplot(1,2,1)
    plot(X,'ro'); %Plot
    subplot(1,2,2)
    plot(X2,'ro'); %Plot
    %drawnow %Als je het nu direct geplot wil
    
    F(j) = getframe(gcf); %Save frame
end


%fig = figure;
%movie(fig,F,2)

v = VideoWriter('testVideo.avi'); %Maak een video-file
open(v)
writeVideo(v,F) %Sla de frames op in de video
close(v)