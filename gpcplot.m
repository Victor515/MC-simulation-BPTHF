function gpcplot(weight,color)%����A�����趨
global t;
global threedplot;
global t_record;
A=2; %���Ƿ�����ת����ͼ����һ�������������޸ģ�����������ͬ����ͬ
m=23-[4:log(1.025315):19]; %��������ϵ
for i=1:1:size(weight,2)
    Mp=24-log(weight(i)); %ת�����������
    wplot(i,:)=1./((A*sqrt(pi))).*exp(-(1/(A^2))*(m-Mp).^2) .*weight(i);

end
wplot = sum(wplot,1);
% mplot=(wplot'*ones(i,1))'; %��wplot��ĳһ��ʱ�̵�gpcͼ��ļ�¼����ת��
% threedplot=[threedplot;mplot];  %�����ת�ú�ľ���ϲ�������¼��һ����ά����
% t_record=[t_record,t];  %ʱ����
% fliplr(t_record);
% [coupleable,couple_time]=coupling_condition();
% if coupleable==0  %��couple�Ѿ��������ͼ
%    [x,y]=meshgrid(m,t_record); %3Dͼ�����Ὠ��
%     mesh(x,y,threedplot);%����
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