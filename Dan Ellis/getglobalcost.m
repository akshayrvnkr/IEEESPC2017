function [pos,C,B] = getglobalcost(x,bpm)
T0=bpm;
L=round(length(x)/bpm);
thresh=T0/5;

C=Inf+zeros(L,length(x));
B=Inf+zeros(L,length(x));

ll(1)=1;
ul(1)=length(x)/L+thresh;
for i=2:L
    ll(i)=round(max(1,ll(i-1)+T0));
    ul(i)=length(x);
end;
for K=1:ul(1)
    n=1:K;
    t=repmat(x(K),length(n),1)';
    t2=n==K;
    C(1,K)=sum(sum((x(n)-t.*t2).^2));
    if(L==1)
        n=1:length(x);
        t=repmat(x(K),length(n),1);
        t2=n==K';
        C(1,K)=sum(sum((x(n)-t.*t2).^2));
    end;
    B(1,K)=0;
end;
for l=2:L
    for K=ll(l):ul(l)
        for m=round(max(1,K-T0-thresh):min(length(x),K-T0+thresh))
            if(m~=K)
                locost=C(l-1,m)+getlocalcost(x,m,K,K);
                if(l==L)
                    locost=C(l-1,m)+getlocalcost(x,m,K,length(x));
                end;
                if(locost<C(l,K))
                    C(l,K)=locost;
                    B(l,K)=m;
                end;
            end;
        end;
    end;
end;
[~,pos(1)]=min(C(end,:));
ctr=2;
for i=L:-1:2
    pos(ctr)=B(i,pos(ctr-1));
    ctr=ctr+1;
end
pos=pos(end:-1:1);
end
function locost = getlocalcost(x,gciloc1,gciloc2,End)
n=gciloc1+1:End;
t1=repmat(x(gciloc2),length(n),1)';
t2=repmat((n==gciloc2),1,1);
locost=sum(sum((x(n)-t1.*t2).^2));
end