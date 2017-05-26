function gpcplot(weight,color)%常数A自行设定
global t;
global threedplot;
global t_record;
A=2; %这是分子量转化成图谱中一个参数，可以修改，随着仪器不同而不同
m=23-[4:log(1.025315):19]; %建立坐标系
for i=1:1:size(weight,2)
    Mp=24-log(weight(i)); %转化成流出体积
    wplot(i,:)=1./((A*sqrt(pi))).*exp(-(1/(A^2))*(m-Mp).^2) .*weight(i);

end
wplot = sum(wplot,1);
% mplot=(wplot'*ones(i,1))'; %将wplot即某一个时刻的gpc图像的记录向量转置
% threedplot=[threedplot;mplot];  %将这个转置后的矩阵合并起来记录成一个二维矩阵
% t_record=[t_record,t];  %时间轴
% fliplr(t_record);
% [coupleable,couple_time]=coupling_condition();
% if coupleable==0  %当couple已经结束后绘图
%    [x,y]=meshgrid(m,t_record); %3D图坐标轴建立
%     mesh(x,y,threedplot);%绘制
% end
switch color
   case 1,
     plot(m,wplot,'-r','LineWidth',0.75);
   case 2,
     plot(m,wplot,'--b');
   case 3,
     plot(m,wplot,':k');
   otherwise
       fprintf('bug');
end