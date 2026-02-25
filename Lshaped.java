package stochasticprograms;

import ilog.cplex.*;
import ilog.concert.*;
public class Lshaped {
	

	//随机规划 L-shaped method
	
	/*
	 * min 10x_1+15x_2+ E[q_1 y_1+q_2 y_2]
	 * x_1+x_2<=12 6y_1+10y_2<=60x_1 8y_1+5y_2<=80x_2
	 * y_1<=d_1 y_2<=d_2 x_1>=4 x_2>=2 y_1,y_2>=0
	 * 
	 * 两种情景  (d_1,d_2,q_1,q_2)
	 * 			0.4, (50,10,-2.4,-2.8)
	 * 			0.6, (30,30,-2.8,-3.2)
	 */
	
	//主问题
	IloCplex mp;
	IloNumVar[] x;
	IloNumVar theta;
	double[] x_value;
	double theta_value;
	public void buildMP() throws IloException
	{
		mp=new IloCplex();
		mp.setOut(null);
		//变量
		x=mp.numVarArray(2, 0, Double.MAX_VALUE);
		theta=mp.numVar(-100000,0);
		//约束
		double[] c1= {1.0,1.0};
		mp.addLe(mp.scalProd(c1, x),12);//x_1+x_2<=12
		mp.addGe(x[0],4);
		mp.addGe(x[1],2);
		//目标
		double[] c= {10.0,15.0};
		IloNumExpr obj=mp.sum(mp.scalProd(c, x),theta);
		mp.addMinimize(obj);
	}
	
	
	//子问题
	//[d_1,d_2,q_1,q_2]
	double[][] dk= {{50,10},{30,30}};
	double[][] qk= {{-2.4,-2.8},{-2.8,-3.2}};
	
	IloCplex sp;
	IloNumVar[] y;
	IloRange[] r;
	double[] p= {0.4,0.6};//概率
	int ncons=4;//约束数
	int scen=2;//2情景
	double[][] h= {{0,0,50,10},{0,0,30,30}};
	
	double[][] T= {{-60,0,0,0},{0,-80,0,0}};
	public void  buildSP(double[] d,double[] q) throws IloException
	{
		sp=new IloCplex();
		sp.setOut(null);
		y=sp.numVarArray(2, 0, Double.MAX_VALUE);
		
		r=new IloRange[4];
		double[] c1= {6.0,10.0};
		double[] c2= {8.0,5.0};
		r[0]=sp.addLe(sp.scalProd(c1,y),60*x_value[0]);
		r[1]=sp.addLe(sp.scalProd(c2,y),80*x_value[1]);
		r[2]=sp.addLe(y[0],d[0]);
		r[3]=sp.addLe(y[1],d[1]);
		
		IloNumExpr obj=sp.scalProd(q,y);
		sp.addMinimize(obj);
		
	}
	double[][] u;
	public void solve() throws IloException
	{
		
		x_value=new double[2];
		u=new double[scen][ncons];//2情景 4约束
		buildMP();
		while(true)
		{
			
			if(mp.solve())
			{
				x_value[0]=mp.getValue(x[0]);
				x_value[1]=mp.getValue(x[1]);
				theta_value=mp.getValue(theta);
			}else
			{
				System.out.println("主问题不可行");
				break;
			}
			
			//对每个情景求解一次
			double e=0;
			double[][] piT=new double[2][2];//对应约束
			double[] E=new double[2];
			double Ex=0;
			for(int i=0;i<scen;i++)
			{
				buildSP(dk[i],qk[i]);
				double pih=0;
				
				if(sp.solve())
				{
					//获取对偶变量u 
					for(int j=0;j<ncons;j++)
					{
						u[i][j]=sp.getDual(r[j]);
						pih+=u[i][j]*h[i][j];
						piT[i][0]+=u[i][j]*T[0][j];
						piT[i][1]+=u[i][j]*T[1][j];
						
					}
					
					
				}else
				{
					System.out.println("子问题求解出错");
				}
				e+=p[i]*pih;
				E[0]+=p[i]*piT[i][0];
				E[1]+=p[i]*piT[i][1];
				
				sp.clearModel();
			}
			Ex=E[0]*x_value[0]+E[1]*x_value[1];
			System.out.println("e="+e+",Ex="+Ex);
			System.out.println("theta="+theta_value);
			double w=e-Ex;
			if(theta_value>=w)
			{
				System.out.println("此时已求得最优");
				System.out.println("x[0]="+x_value[0]+",x[1]="+x_value[1]);
				break;
			}else
			{
				//添加最优割到主问题
				IloRange cut=mp.ge(mp.sum(theta,mp.scalProd(E, x)),e);
				mp.addCut(cut);
			}
			
			
		}
		
	}
	public static void main(String[] args) throws IloException {
		Lshaped ls=new Lshaped();
		ls.solve();

	}

}
