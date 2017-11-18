#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <pthread.h>
#define rep(i,x,y) for(int i=x;i<y;i++)
#define N 4000000
#define B2M 8000000
#define NUM_THREADS     1
using namespace std;
typedef long long LL;

struct matrix{
    LL x;int I;
}G[B2M],g[B2M];
struct matrix2{
    LL x;int I,J;
}G2[B2M];
struct equation{
    int q1,q2,q3,q4;
    LL pat;
}F[B2M],F1[B2M];
struct t_data{
    int begin;
};
//int B2=40,B=14,Round=1,top=59-Round*B,bot=60-(Round+1)*B;
int B2=40,B=14,Round=1,top=51,bot=38;
LL bi[64];
LL m=0;
LL M=0;
const int max_E=16385;
int E[max_E],Sum[max_E],z[N];
int v=8;

void init_bi()
{
    bi[0]=1;
    rep(i,1,64)
        bi[i]=bi[i-1]*2;
}

LL miss(LL x,int top,int bot)
{
    return bi[bot]*(x/bi[top+1])+x%bi[bot];
}
LL got(LL x,int top,int bot)
{
    return (x%bi[top+1])/bi[bot];
}
LL up(LL x,LL top, LL bot)
{
    return miss(x,top,bot)+got(x,top,bot)*bi[60-B];
}
bool cmpB2(matrix A,matrix B)
{
    return A.x%bi[60-B2-v]>B.x%bi[60-B2-v];
}
bool cmpB(matrix2 A,matrix2 C)
{
    return A.x%bi[60-B]>C.x%bi[60-B];
}
bool cmpM(equation A,equation C)
{
    if(A.pat<C.pat)
        return true;
    else if(A.pat>C.pat)
        return false;
    else if(A.q1<C.q1)
        return true;
    else if(A.q1>C.q1)
        return false;
    else if(A.q2<C.q2)
        return true;
    else if(A.q2>C.q2)
        return false;
    else if(A.q3<C.q3)
        return true;
    else if(A.q3>C.q3)
        return false;
    else if(A.q4<C.q4)
        return true;
    else return false;
}
bool eql(equation A,equation C)
{
    if(A.pat==C.pat && A.q1==C.q1 && A.q2==C.q2 && A.q3==C.q3 && A.q4==C.q4)
        return true;
    else
        return false;
}
void bin(int x)
{
    int binx[60]={0};
    int top=0;
    cout<<"x="<<x<<endl;
    while(x>0)
    {
        if(x%2==0) binx[top++]=0;
        else       binx[top++]=1;
        x/=2;
    }
    for(int i=B-1;i>=0;i--)
        cout<<binx[i];
    cout<<endl;
}
void init_G(){
    G[0].x=pow(2,59);

    for(int i=1;i<60;i++)
    {
        G[i].x=G[i-1].x/2;
        G[i].I=i;
        
    }
    rep(i,0,60)
    {
        
        G[i].x=up(G[i].x,top,bot);
        
        LL res=G[i].x;
        /*
        cout<<59-i<<":\t";
        while(res)
        {
            cout<<(res&1);
            res/=2;
        }
        cout<<endl;*/
    }
    for(int i=60;i<N;i++){
        G[i].x=G[i-3].x^G[i-9].x^G[i-16].x^G[i-35].x^G[i-37].x^G[i-47].x^G[i-56].x^G[i-60].x;
        G[i].I=i;
    }

}
void init_g(){
    g[0].x=pow(2,59);

    for(int i=1;i<60;i++)
    {
        g[i].x=g[i-1].x/2;
        g[i].I=i;
    }

    for(int i=60;i<N;i++){
       g[i].x=g[i-3].x^g[i-9].x^g[i-16].x^g[i-35].x^g[i-37].x^g[i-47].x^g[i-56].x^g[i-60].x;
        g[i].I=i;
    }

}
void init_z(char data[]){
    FILE* f=fopen(data,"r");
    rep(i,0,N)
        fscanf(f,"%d",&z[i]);
    fclose(f);
}
void a1_step1(int vary)
{
    m=0;
    //rep(i,0,100)
        //cout<<G[i].x<<endl;
      //  bin(G[i].x);
    sort(G+60,G+N,cmpB2);
    rep(i,60,N)
    {
        int j=i;
        while(((G[j].x^G[j+1].x)%bi[60-B2-v])==0 && j<N) j++;
        if(i==j)    continue;
        
                        
        rep(k1,i,j)
            rep(k2,k1+1,j+1)
            {


                LL X=G[k1].x^G[k2].x;
                if(got(X,19,20-v)==vary){
                    G2[m].x=X;
                    G2[m].I=G[k1].I;G2[m].J=G[k2].I;
                    m++;
                    
                }
            }
        i=j;
    }
    cout<<"m="<<m<<endl;
    sort(G2,G2+m,cmpB);
    
    rep(i,60,m)
    {
        int j=i;
        while(G2[j].x%bi[60-B]==G2[j+1].x%bi[60-B] && j<m) j++;
        if(i==j)    continue;

        rep(k1,i,j)
            rep(k2,k1+1,j+1)
            {
                int i1=G2[k1].I,j1=G2[k1].J,i2=G2[k2].I,j2=G2[k2].J;

                if(i1==i2 || i1==j2 || j1==i2 || j1==j2) continue;


                LL e0=G2[k1].x^G2[k2].x;

                e0/=bi[60-B];
                 E[e0]+=1;
                Sum[e0]+=z[i1]^z[i2]^z[j1]^z[j2];
                M++;
            }
        i=j;
    }

}
int num_in_bin(int x)
{
    int ones=0;
    while(x>0)
    {
        if(x%2==1) ones^=1;
        x/=2;
    }
    return ones;
}

void a2_step1()
{

    cout<<"old M is "<<M<<endl;
    sort(F,F+M,cmpM);
    int num=0;
    rep(i,0,M)
    {
        int j=i;
        while(eql(F[j],F[j+1])) j++;
        F1[num].pat=F[i].pat;      F1[num].q1=F[i].q1;      F1[num].q2=F[i].q2;     F1[num].q3=F[i].q3;     F1[num].q4=F[i].q4;
        num++;
        i=j;
    }

    cout<<"new M is "<<num<<endl;

    rep(i,0,num)
    {
        E[F1[i].pat]+=1;
        Sum[F1[i].pat]+=z[F1[i].q1]^z[F1[i].q2]^z[F1[i].q3]^z[F1[i].q4];
    }
    return ;
}
void a2_step2()
{

    int u_point,max_point=0,best_guess=-1;
    rep(u,0,bi[B])
    {
        u_point=0;
        rep(e,0,bi[B])
        {
            int Left = num_in_bin(e&u);
            if(Left==1) u_point+=Sum[e];
            else        u_point+=E[e]-Sum[e];

        }
        if(u_point>max_point){
            max_point=u_point;
            best_guess=u;

        }
    }
    bin(best_guess);
    printf("\nM=%lld,m_p=%d,per=%.3lf\n",M,max_point,(double)max_point/(double)M);
}
int main(int argv,char* argc[])
{
    top=(argc[1][0]-'0')*10+argc[1][1]-'0';
    bot=(argc[2][0]-'0')*10+argc[2][1]-'0';
    char data[20]="data-09.txt";
	data[5]=argc[3][0];
	data[6]=argc[3][1];		    
    if(strlen(argc[4])==2)
   	 v=(argc[4][0]-'0')*10+argc[4][1]-'0';
    else
	v=argc[4][0]-'0';   
    cout<<top<<endl;
    cout<<bot<<endl;
    B=top-bot+1;
    init_bi();
    init_G();
    //init_g();
    init_z(data);
    memset(E,0,sizeof(E));

    rep(i,0,bi[v]-1)     a1_step1(i);
    rep(i,0,59-top) cout<<"_";
    a2_step2();
    rep(i,0,bot) 	cout<<"_";
    return 0;
}
