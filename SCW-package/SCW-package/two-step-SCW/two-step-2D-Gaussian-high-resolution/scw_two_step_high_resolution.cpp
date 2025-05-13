///remove all zero lines in theta matrix; 
/// and add 0.2tra + 0.8ter, in this version change this, using the count ratio with inter/intra
// this will add expend list content
//expended criteria is followed scl
//
//
//
//using 10m - 1m, same criteria for high resolution, see resutls
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <sys/time.h>
#include <thread>
#include <algorithm>
#include <unordered_map>
#include <string.h>

using namespace std;

int size_c;
int size_c2;
int np;
double res;
double res1;
double res2;
double d0 = 8.0;
int size_b;
int size_b2;//for expand_list 
double c_ratio;//contacts raiton intra/inter;
vector <int> len_chr;
vector <int> len_chr2;

vector <double> cube_cost;
vector <vector <string>> cube_io;
vector <vector <int>> cube_id;
vector <vector <double>> slist;
vector <vector <double>> slist1;
vector <vector <double>> slist2;
vector <vector <vector <double>>> cube_list;

vector <string> ordered_chr;
vector <string> list_chr;
vector <int> list_chr_id;
vector <string> list_chr2;
vector <int> list_chr_id2;
vector <vector <string>> theta_list;

unordered_map <string, double> thEta;
unordered_map <string, double> thEta2;
unordered_map <string, double> HiC;
unordered_map <string, double> HiC2;
//unordered_map <string, int> map_key;
unordered_map <string, int> map_del;
unordered_map <string, int> map_del2;
unordered_map <int, string> key_map;
unordered_map <int, string> key_map2;
unordered_map <int, string> del_map;
unordered_map <int, string> del_map2;
unordered_map <string, int> scaled_map;

vector <vector <string>> readMatrixFile(char* file){
	vector <vector <string>> matrix;
	ifstream myfile;
	myfile.open(file);
	string line;
	string temp;
	while(getline(myfile, line)){
		vector <string> row;
		istringstream istr(line);
		while (istr >> temp){
			string tempf = temp.c_str();
			row.push_back(tempf);
		}
		matrix.push_back(row);	
	}
	return matrix;
}




unordered_map <string,double >readHic(vector<vector <string>> mat, double reso){
	int count1 = 0;
	int count2 = 0;
	unordered_map <string, double> hMat;
	for(int i = 0; i < mat.size(); i++){
		int s1 = int(stoi(mat[i][1])/reso);
		int s2 = int(stoi(mat[i][3])/reso);
		string name = mat[i][0]+" "+to_string(s1)+" " + mat[i][2] + " " + to_string(s2);
		string name2 = mat[i][2]+" "+to_string(s2)+" " + mat[i][0] + " " + to_string(s1);
		//cout<<name<<endl;
		//cout<<name2<<endl;
		if(mat[i][0] == mat[i][2]){
			count1++;
		}else{
			count2++;
		}
		hMat[name] = 1;
		hMat[name2] = 1;
	}
	//cout<<"Inter contacts: "<<count2<<"; Intra contacts: "<<count1<<endl;
	c_ratio = int(count1/count2);
	c_ratio = 1/(c_ratio + 1);
	//cout<<c_ratio<<endl;
	return hMat;
}


void calTheta(vector<vector <int>> ones, int m, int n,int count, vector <int> clen){
	//vector <vector <int>> ones;
	//cout<<m<<" "<<n<<" "<<ones.size()<<endl;

	for(int i = 0; i < clen[m]; i++){
                for(int j = 0; j < clen[n]; j++){
                        string name = ordered_chr[m]+" "+to_string(i)+" "+ordered_chr[n]+" "+to_string(j);
			string name2= ordered_chr[n]+" "+to_string(j)+" "+ordered_chr[m]+" "+to_string(i);
			//cout<<name;
			double temp1 = 0;
			int is0 = 0;
                        for(int k = 0; k < ones.size(); k++){
				if((ones[k][0] == i) && (ones[k][1] == j)){
					is0 = 1;
				}
			}
			if(is0 == 0){
				for(int k = 0; k < ones.size(); k++){
					double i1 = ones[k][0];
					double j1 = ones[k][1];
					double temp2 = (i1-i)*(i1-i)+(j1-j)*(j1-j);
					temp2 = temp2/2;// this two is mu2,
					//cout<<temp2<<" ";
					temp2 = exp(-temp2);
					//cout<<temp2<<" ";
					temp1 += temp2;
					//cout<<" "<<temp2;	
				}
				//cout<<temp1<<endl;	
				if(temp1 >= 1){
					temp1 = 1;
				}
				if(temp1 > 0.7){
					theta_list[count].push_back(name+" "+to_string(temp1));
					theta_list[count].push_back(name2+" "+to_string(temp1));
					//thEta[name] = temp1;/// could test to change this 0.7// this may met pthread problem when size is small, 
					//thEta[name2]= temp1;
				}
				//cout<<endl;
				//cout<<name<<" "<<temp1<<endl;
			}else{
				temp1 = 1;
				theta_list[count].push_back(name+" "+to_string(temp1));
				theta_list[count].push_back(name2+" "+to_string(temp1));
					//thEta[name] = 1;
					//thEta[name2] = 1;				
			}
                        
                }
        }
	//cout<<m<<" "<<n<<endl;
	
}


unordered_map <string, double> readTheta(unordered_map <string, double> hMat, vector <int> clen){
	const size_t THREAD_COUNT = np;//number of thread
	unordered_map <string, double> tMat;
	theta_list.clear();
	vector <thread> threadPool;
	for(int i = 0; i < ordered_chr.size(); i++){
		for(int j = i; j < ordered_chr.size(); j++){
			vector <string> tempa;
			theta_list.push_back(tempa);
	
		}
	}

	int count = 0;
	for(int i = 0; i < ordered_chr.size(); i++){
		for(int j = i; j < ordered_chr.size(); j++){

			vector <vector <int>> ones;
     
	        	for(int i1 = 0; i1 < clen[i]; i1++){
         	        	for(int j1 = 0; j1 < clen[j]; j1++){
                	        	string name = ordered_chr[i]+" "+to_string(i1)+" "+ordered_chr[j]+" "+to_string(j1);
                        		//cout<<name<<" "<<HiC[name]<<endl;
					if(hMat.count(name)){
						//cout<<name<<" "<<i1<<" "<<j1<<" "<<HiC[name]<<endl;
                                		//thEta[name] = 1;
                                		vector <int> row;
                                		row.push_back(i1);
                                		row.push_back(j1);
                                		ones.push_back(row);
                        		}
                		}
        		}
			//if(ones.size() != 0){
			threadPool.emplace_back([ones,i,j,count,clen](){
				calTheta(ones,i,j,count,clen);				
			}); 
			//}
			count++;
		}
	}
	for(auto& t: threadPool){
		t.join();
	}

        for(int i = 0; i < theta_list.size(); i++){
                for(int j = 0; j < theta_list[i].size(); j++){
			//cout<<theta_list[i][j]<<endl;
			
               		vector <string> row;
			string temp;
			istringstream istr(theta_list[i][j]);
			while (istr >> temp){
				string tempf = temp.c_str();
				row.push_back(tempf);
			}
			string name = row[0]+" "+row[1]+" "+row[2]+" "+row[3];
			tMat[name] = stod(row[4]);
                }
        }

	return tMat;
		
	
}



/*void see_list(){
	for(int i = 0; i < len_chr.size(); i++){
		int count = 0;
		for(int j = 0; j < len_chr[i]; j++){
			string name = ordered_chr[i] + " " + to_string(j);
			if(map_del.count(name)){
				count++;
			}
		}
		cout<<ordered_chr[i]<<" "<<ct<<endl;
	}
}*/















vector <vector <double>> find_ninit(vector <vector <double>> lis, int sizec){
	//cout<<lis.size()<<endl;
	vector <vector <double>> possible_pos;
	double x_bottom_bound =((double)lis[lis.size()-1][0])-sqrt(10)-1;
	double x_upper_bound = ((double)lis[lis.size()-1][0])+sqrt(10)+1;
	for(int x = (int)x_bottom_bound; x <= (int)x_upper_bound;x++){
		if(x<0||x>=sizec){
			continue;
		}
		double y_bottom_bound =((double)lis[lis.size()-1][1])-sqrt(10)-1;
        	double y_upper_bound = ((double)lis[lis.size()-1][1])+sqrt(10)+1;
		for(int y = (int)y_bottom_bound; y <= (int)y_upper_bound;y++){
			if(y<0||y>=sizec){
				continue;
			}
			double z_bottom_bound =((double)lis[lis.size()-1][2])-sqrt(10)-1;
                	double z_upper_bound = ((double)lis[lis.size()-1][2])+sqrt(10)+1;
			for(int z = (int)z_bottom_bound; z <= (int)z_upper_bound;z++){
				if(z<0||z>=sizec){
					continue;
				}
				//cout<<x<<" "<<y<<" "<<z<<endl;
				vector <double> lll1;
				lll1.push_back(x);
				lll1.push_back(y);
				lll1.push_back(z);
				possible_pos.push_back(lll1);
				
			}
		}
	}
	vector <vector <double>> possible_pos2;
	bool random_find = false;
	int count = 0;
	while(random_find == false){
		count++;
		if(count==possible_pos.size()){
			break;
		}
		timeval t1;
		gettimeofday(&t1,NULL);
		srand(t1.tv_usec * t1.tv_sec);
		int rand1 = rand() % possible_pos.size();
		//cout<<rand1<<endl;
		int x_p = possible_pos[rand1][0];
		int y_p = possible_pos[rand1][1];
		int z_p = possible_pos[rand1][2];
		//cout<<x_p<<
		bool removed = false;
		for(int i = 0; i < slist.size();i++){
			double x_e = slist[i][0];
			double y_e = slist[i][1];
			double z_e = slist[i][2];
			double dist = sqrt((x_p-x_e)*(x_p-x_e)+(y_p-y_e)*(y_p-y_e)+(z_p-z_e)*(z_p-z_e));
			if(i < slist.size()-1){
				if(dist < 2 || dist == sqrt(8)){
					removed = true;
					break;
				}
			}
			if(i==slist.size()-1){
				if(dist < 2 || dist == sqrt(8)||dist>sqrt(10)){
					removed = true;
					break;
				}
			}
		}
		if(removed == false){
			possible_pos2.push_back(possible_pos[rand1]);
		}
	}
	if(possible_pos2.size()<0){
		exit(1);
	}
	return possible_pos2;
}




bool rand_init(int len1, int len2){
	slist.clear();
        for(int i = 0; i < len1; i++){
                if(i == 0){
                        timeval t1;
                        gettimeofday(&t1, NULL);
                        srand(t1.tv_usec * t1.tv_sec);

                        int rand1 = rand() % len2;
                        int rand2 = rand() % len2;
                        int rand3 = rand() % len2;
                        vector <double> sl;
                        sl.push_back(rand1);
                        sl.push_back(rand2);
                        sl.push_back(rand3);
                        slist.push_back(sl);
                }else{
                        vector <vector <double>> candidates = find_ninit(slist, len2);
                        if(candidates.size() == 0){
                                cout<<"The size of candidate nodes are 0! initialization failed! There are no enough space. Erase everything and restart"<<endl;
                                slist.clear();
                                return false;
                        }

                        slist.push_back(candidates[0]);
                }
	}
        return true;
}



/*void mid(){
	double x, y, z;
	for(int i = 0; i < size_b; i++){
		x += slist[i][0];
		y += slist[i][1];
		z += slist[i][1];
	}
	x /= size_b;
	y /= size_b;
	z /= size_b;
	double dx = size_c/2 - x;
	double dy = size_c/2 - y;
	double dz = size_c/2 - z;
	for(int i = 0; i < size_b; i++){
		slist[i][0] += dx;
		slist[i][1] += dy;
		slist[i][2] += dz;
	}

}*/


void gen_list(int len1, vector <int> id, vector <string> io, vector <vector <double>> cords){
	cube_id.clear();
	cube_io.clear();
	cube_list.clear();
	cube_cost.clear();	
	int len2 = len1/np + 1;
	int cube_num = int(len1/len2) + 1;
	for(int i = 0; i < cube_num; i++){
		vector<vector <double>> l1;
		cube_list.push_back(l1);
		vector<string> l2;
		cube_io.push_back(l2);
		vector<int> l3;
		cube_id.push_back(l3);
	}	
	for(int i = 0; i < cube_num; i++){
		int len3 = len2; //end of the nums;
		if(i == cube_num - 1){
			len3 = len1 - (cube_num-1)*len2;
		}
		vector <vector <double>> clist;
		vector <int> cid;
		vector <string> cio; 
		for(int j = 0; j < len3; j++){
			int k = i*len2 + j;
			//cout<<i<<" "<<k<<"; Chr"<<list_chr[k]<<" : "<<list_chr_id[k]<<endl;
			cid.push_back(id[k]);
			cio.push_back(io[k]);
			clist.push_back(cords[k]);
		}
		cube_id[i] = cid;
		cube_io[i] = cio;
		cube_list[i] = clist;
		cube_cost.push_back(0);
	}
	//cout<<len2<<endl;
}


void cal_cube(int num, int m, int n, bool cal){
	//cout<<num<<" "<<m<<" "<<n<<endl;
	//double cost1 =0.0;
	//double cost1_1=0.0;
	//double cost1_2=0.0;
	//double cost1_3=0.0;
	double cost2 = 0.0;
        double t2 = 0;
        t2 = pow(0.7,1.0/3);
        t2 = d0/t2;
	double xi = cube_list[m][n][0];
	double yi = cube_list[m][n][1];
	double zi = cube_list[m][n][2];
	for(int i = 0; i < cube_list[num].size();i++){
	    if(m != num || n != i){
		double xj = cube_list[num][i][0];
		double yj = cube_list[num][i][1];
		double zj = cube_list[num][i][2];
		double cost1 = 0.0;
		double cost1_1 = 0.0;
		double cost1_2 = 0.0;
		double cost1_3 = 0.0;
		double dist = sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj));
		string name = cube_io[m][n] + " " + to_string(cube_id[m][n]) + " " + cube_io[num][i] + " "+to_string(cube_id[num][i]);
		//cout<<name<<endl;
		double tva;
		int yes = 0;
		if(cal == false){
			if(thEta.count(name)){
				yes = 1;
				tva = thEta[name];
			}
		}else{
			//tva = thEta2.count(name);
			if(thEta2.count(name)){
				yes = 1;
				tva = thEta2[name];
			}
		}
		if(yes == 1){
			if(tva == 1){
				//if(dist < d0){	
					cost1_1 = cost1_1 + (((dist - d0) * (dist - d0))/(d0*d0));//this 2 is target values,
				//}
			}
			if(tva >= 0.7 && tva < 1){
				double t1 = tva;
				t1 = pow(t1, 1.0/3);
				t1 = d0/t1;//edit 2 as para
				//if(dist < t1){
					cost1_2 = cost1_2 + 1*(1 - exp(-((dist-t1)*(dist-t1))/20.0));//fist 1 is beta, 20 is mu1;
				//}
			}
			/*if(thEta[name] < 0.7){
				double t2 = 0;
				t2 = pow(0.7,1.0/3);
				t2 = 2/t2;//this 2 could be changed
				cost1_3 += 1*(1-1/(1+exp(-(dist-(t2-1)))/0.1));
			}*/
		}else{
			//if(thEta[name] < 0.7){
                        //double t2 = 0;
                        //t2 = pow(0.7,1.0/3);
                        //t2 = 8.0/t2;//this 2 could be changed
                        //if(dist < t2){
			cost1_3 = cost1_3 + 1*(1-1/(1+exp(-(dist-(t2-1.0)))/0.1));
                        //}
			//}*/
			
		}
		cost1 = cost1_1 + cost1_2 + cost1_3;
                if(cube_id[m][n] == cube_id[num][i]){
                        cost1 = (c_ratio) * cost1;
                }else{
                        cost1 = (1-c_ratio) * cost1;
                }
                cost2 = cost2 + cost1;
	    }		
	}
	//cost1 = cost1_1 + cost1_2 + cost1_3;
	cube_cost[num] = cost2;
}
	

void call_thread(int m, int n, bool cal){
	const size_t THREAD_COUNT = cube_cost.size();
	vector <thread> threadPool;
	for(int i = 0; i < cube_cost.size(); i++){
		threadPool.emplace_back([i,m,n,cal](){
			cal_cube(i,m,n,cal);
		});
	}
	for(auto& t: threadPool){
                t.join();
        }
}


double calculate_loss(int m, int n, bool cal){
	double cost = 0.0;
	for(int i = 0; i < cube_cost.size(); i++){
		cube_cost[i] = 0;
	}
	call_thread(m, n, cal);
	for(int i = 0; i < cube_cost.size(); i++){
		cost = cost + cube_cost[i];
	}
	return cost;
}


bool simulate_annealing(double temperature, unordered_map <string, int> index, bool cal){
	timeval t1;
	gettimeofday(&t1, NULL);
	srand(t1.tv_usec * t1.tv_sec);
	int rand1 = rand()%cube_list.size();
	int rand2 = rand()%cube_list[rand1].size();
	double cost1 = calculate_loss(rand1, rand2, cal);
	//cout<<cost1<<endl;	
	double oldx = cube_list[rand1][rand2][0];
	double oldy = cube_list[rand1][rand2][1];
	double oldz = cube_list[rand1][rand2][2];
	//cout<<"Old: "<<oldx<<" "<<oldy<<" "<<oldz<<" Cost1: "<<cost1<<endl;
	string name = cube_io[rand1][rand2] + " " + to_string(cube_id[rand1][rand2]);
	int cur = index[name];
	vector <vector <double>> cand1;
        double mod = 1;
        for(double i = oldx-mod;i<= oldx+mod; i=i+mod){
                for(double j = oldy -mod; j <= oldy+mod; j=j+mod){
                        for(double k = oldz -mod; k <= oldz+mod;k=k+mod){
				if(i == oldx && j == oldy && k == oldz){
					continue;
				}
				else{
					if(cur == 0){
						double xmax = slist[1][0];
						double ymax = slist[1][1];
						double zmax = slist[1][2];
						double dist2 = sqrt((xmax - i)*(xmax - i)+(ymax-j)*(ymax-j)+(zmax-k)*(zmax-k));
						if(dist2 > 8){
							continue;
						}
					}else if(cur == slist.size() - 1){
						double xmin = slist[cur-1][0];
						double ymin = slist[cur-1][1];
						double zmin = slist[cur-1][2];
						double dist1 = sqrt((i-xmin)*(i-xmin)+(j-ymin)*(j-ymin)+(k-zmin)*(k-zmin));
						if(dist1 >8){
							continue;
						}
					}else{
						double xmax = slist[cur+1][0];
                                                double ymax = slist[cur+1][1];
                                                double zmax = slist[cur+1][2];
                                                double dist2 = sqrt((xmax - i)*(xmax - i)+(ymax-j)*(ymax-j)+(zmax-k)*(zmax-k));
                                                double xmin = slist[cur-1][0];
                                                double ymin = slist[cur-1][1];
                                                double zmin = slist[cur-1][2];
                                                double dist1 = sqrt((i-xmin)*(i-xmin)+(j-ymin)*(j-ymin)+(k-zmin)*(k-zmin));
						if(dist1 > 8 || dist2 > 8){
							continue;
						}
	
					}
						
					
                                	vector <double> xyz_c;
                                	xyz_c.push_back(i);
                                	xyz_c.push_back(j);
                                	xyz_c.push_back(k);
                                	cand1.push_back(xyz_c);
				}
				
                        }
                }
        }
	if(cand1.size() == 0){
		return false;
	}
        int rand3 = rand()%cand1.size();
	double newx = cand1[rand3][0];
	double newy = cand1[rand3][1];
	double newz = cand1[rand3][2];
	//cout<<"New: "<<newx<<" "<<newy<<" "<<newz<<endl;

        cube_list[rand1][rand2][0] = newx;
        cube_list[rand1][rand2][1] = newy;
        cube_list[rand1][rand2][2] = newz;
	slist[cur][0] = newx;
	slist[cur][1] = newy;
	slist[cur][2] = newz;

	double cost2 = calculate_loss(rand1,rand2, cal);
	//cout<<"rand cands has "<<cand1.size()<<" Picked "<<rand3<<" th."<<endl;
	//cout<<"New: "<<newx<<" "<<newy<<" "<<newz<<": cost2 is "<<cost2<<endl;

	double Z = cost2 - cost1;
	if(Z <= 0){
		//cout<<" Z <= 0 ture"<<endl<<endl;
		//cout<<"Loss "<<cost2<<endl;
		return true;
	}else{
		double rand4 = 0;
                timeval t2;
                gettimeofday(&t2,NULL);
                srand(t2.tv_usec * t2.tv_sec);
		rand4 = rand()%(999999+1)/(float)(999999+1);
		//cout<<"random a number "<<rand4<<" Temerature is "<<temperature<<endl;
		if(rand4 < exp(-Z/temperature)){
		
			//cout<<"random number is "<<rand4<<" temperature is "<<temperature;
			//cout<<" check with "<<exp(-Z/temperature)<<" ; less true"<<endl<<endl;
			//cout<<"Loss "<<cost2<<endl;
			return true;
		}else{
			cube_list[rand1][rand2][0] = oldx;
		        cube_list[rand1][rand2][1] = oldy;
		        cube_list[rand1][rand2][2] = oldz;
			slist[cur][0] = oldx;
			slist[cur][1] = oldy;
			slist[cur][2] = oldz;
			//cout<<"Loss "<<cost1<<endl;
                        //cout<<"random number is "<<rand4<<" temperate is "<<temperature;
                        //cout<<" check with "<<exp(-Z/temperature)<<" ; large false"<<endl<<endl;
			return false;
		}

	}
}


void expand_list(){
	vector <vector <double>> nlist;
	vector <string> nid;
	double scal = res1/res2;
	double tims = 0.12 * scal;
//	cout<<scal<<endl; 

        int count = 0;
	int count2= 0;

	for(int i = 0; i < slist.size(); i++){
		slist[i][0] *= tims;
		slist[i][1] *= tims;
		slist[i][2] *= tims;
	}

	int countz = 0;
	for(int i = 0; i < slist.size(); i++){
		vector <double> lis;
		lis.push_back(slist[i][0]);
		lis.push_back(slist[i][1]);
		lis.push_back(slist[i][2]);
		double s1, s2, s3;
		//nlist.push_back(lis);
		if(i != slist.size() - 1){
			s1 = (slist[i+1][0] - slist[i][0])/scal;
			s2 = (slist[i+1][1] - slist[i][1])/scal;
			s3 = (slist[i+1][2] - slist[i][2])/scal;
		}else{
			s1 = (slist[i-1][0] - slist[i][0])/scal;
                        s2 = (slist[i-1][1] - slist[i][1])/scal;
                        s3 = (slist[i-1][2] - slist[i][2])/scal;
		}
		int p = 0;
			//countz++;
			//cout<<del_map[i]<<" ";
		istringstream iss(del_map[i]);
		vector <string> tokens;
		string token;
		while(iss>>token){
			tokens.push_back(token);
		}
		//cout<<int(stoi(tokens[1]) * res1)<<endl;
		//cout<<i*scal<<" "<<i*scal + scal<<endl;
		for(int j = i * scal; j < i*scal + scal; j++){
			vector <double> lis2;
			double a1 = slist[i][0] + p*s1;
			double a2 = slist[i][1] + p*s2;
			double a3 = slist[i][2] + p*s3;
			lis2.push_back(a1);
			lis2.push_back(a2);
			lis2.push_back(a3);
			nlist.push_back(lis2);
			string name = tokens[0] + " " + to_string(int(stoi(tokens[1]) * scal + p));
			nid.push_back(name);
			//cout<<tokens[0]<<" "<<int(stoi(tokens[1]) * scal + p)<<endl;	
			p++;
			countz++;
		}
		
	}
	//cout<<countz<<" "<<nlist.size()<<endl;
	/*for(int i = 0; i < size_b2; i++){
		string name3 = del_map2[i];
		if(nid.count(name3)){
	
		}else{
			cout<<name2<<endl;
		}
	}*/
	int countd = 0;
	for(int i = 0; i < nid.size(); i++){
		if(map_del2.count(nid[i])){
			slist2.push_back(nlist[i]);
			//cout<<nid[i]<<" "<<nlist[i][0]<<" "<<nlist[i][1]<<" "<<nlist[i][2]<<endl;
			//countd++;
		}
	}
	//cout<<"Finish extended the list, now has "<<countd<<" beads"<<endl;
	//exit(1);

	
}


void cooling(double temperature, int size, unordered_map<string, int> index, bool cal){
	cout<<"Cooling started."<<endl;
	int failed_times = 0;
	for(int i = 1; i < 1000000; i++){
		cout<<"Temperature: "<<temperature<<endl;
		int counter = 0;
		int accept_counter = 0;
		while(counter < (size * 100)){
			counter++;
			if(simulate_annealing(temperature, index,cal) == true){
				accept_counter++;
			}
			if(accept_counter == (size * 10)){
				failed_times = 0;
				break;
			}
			if(counter == ((size * 100)-1)){
				failed_times++;
			}
			        /*for(int i = 0; i < cube_list.size(); i++){
			                for(int j = 0; j < cube_list[i].size();j++){
                        			cout<<cube_list[i][j][0]<<" "<<cube_list[i][j][1]<<" "<<cube_list[i][j][2]<<" "<<cube_io[i][j]<<" "<<to_string(int(cube_id[i][j]*res2))<<"-"<<to_string(int((cube_id[i][j]+1)*res2))<<endl;
                			}
        			}*/
		}
		//cout<<"Accept num: "<<accept_counter<<endl;
		if(cal == true){
			temperature = temperature * 0.6;
		}else{
			temperature = temperature * 0.9;
		}
		
		if(failed_times == 3){
			break;
		}
		if(cal == true){
			if(temperature < 0.1){
				break;
			}
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////int main //////////////////////
int main(int argc, char* argv[]){
	char* hic_file;
	char* length_file;
	char* out_file;
	for(int i = 1; i < argc - 1; i++){
		if(strcmp(argv[i],"-res") == 0){
			string r (argv[i+1]);
			res = atof(r.c_str());
		}else if(strcmp(argv[i],"-len")==0){
			length_file = argv[i+1];
		}else if(strcmp(argv[i],"-i")==0){
			hic_file = argv[i+1];
		}else if(strcmp(argv[i],"-np")==0){
			string r (argv[i+1]);
			np = atoi(r.c_str());	
		}else if(strcmp(argv[i],"-o")==0){
			out_file = argv[i+1];	
		}
	}
	

	res1 = 1 * 1000000;//basic res
	res2 = res* 1000000;//target res

	vector <vector <string>> clines = readMatrixFile(length_file);
	//cout<<"Chromosome numbers: "<<clines.size()<<endl;
	//int whole_beads = 0;
	int count_beads = 0;// hash{count} = X chr_index
	int count_beads2= 0;
	for(int i = 0; i < clines.size(); i++){
		int num = stoi(clines[i][1]);
		int num2 = stoi(clines[i][1]);
		if(num % int(res1) == 0){
			num = num/(res1);
		}else{
			num  = num/(res1) + 1;
		}
                if(num2 % int(res2) == 0){
                        num2 = num2/(res2);
                }else{  
                        num2  = num2/(res2) + 1;
                }
		for(int j = 0; j < num; j++){
			string name = clines[i][0] + " " + to_string(j);
			key_map[count_beads] = name;
			count_beads++;
		}
		for(int j = 0; j < num2;j++){
			string name = clines[i][0] + " " + to_string(j);
			key_map2[count_beads2] = name;
			count_beads2++;
		}
		size_b += num;
		size_b2 += num2;
		len_chr.push_back(num);
		len_chr2.push_back(num2);
		ordered_chr.push_back(clines[i][0]);	
	}
	cout<<"Beads numbers: "<<"1 Mb "<<size_b<<"; "<<res<<" Mb "<<size_b2<<endl;
	//exit(1);
	vector <vector <string> > hicmat = readMatrixFile(hic_file); 


	//cout<<"Contacts: "<<hicmat.size()<<endl;
	//cout<<"Imputation started!"<<endl;
	HiC = readHic(hicmat, res1);// if met an error of segment core, that is the 0 base exceed problem.
	HiC2 = readHic(hicmat,res2);
	//exit(1);
	int lt = 0;
	int lt2= 0;
	//delete null HiC beads; this will cal tiwice one basic, one target res
	for(int i = 0; i < size_b; i++){
        	int ct = 0; 
		for(int j = 0; j < size_b; j++){
                        string s1 = key_map[i];
                        string s2 = key_map[j];
                        string name = s1 + " " + s2;
                        if(HiC.count(name)){
                        }else{
				ct++;
                        }
                }
		if(ct != size_b){
			del_map[lt] = key_map[i];
			map_del[key_map[i]] = lt;
			vector <string> row;
                        string temp;
                        istringstream istr(key_map[i]);
                        while (istr >> temp){
                                string tempf = temp.c_str();
                                row.push_back(tempf);
                        }
			list_chr.push_back(row[0]);
			list_chr_id.push_back(stoi(row[1]));
			lt++; 
		}
		if(ct == size_b){
		}
        }
	
        for(int i = 0; i < size_b2; i++){
                int ct = 0;
                for(int j = 0; j < size_b2; j++){
                        string s1 = key_map2[i];
                        string s2 = key_map2[j];
                        string name = s1 + " " + s2;
                        if(HiC2.count(name)){
                        }else{
                                ct++;
                        }
                }
                if(ct != size_b2){
                        del_map2[lt2] = key_map2[i];
			//cout<<del_map2[lt2]<<endl;
                        map_del2[key_map2[i]] = lt2;
                        vector <string> row;
                        string temp;
                        istringstream istr(key_map2[i]);
                        while (istr >> temp){
                                string tempf = temp.c_str();
                                row.push_back(tempf);
                        }
                        list_chr2.push_back(row[0]);
                        list_chr_id2.push_back(stoi(row[1]));
                        lt2++;
                }
                if(ct == size_b2){
                }
        }


	size_b = lt;
	size_b2=lt2;
	cout<<lt<<" "<<lt2<<endl;

	//exit(1);
	thEta = readTheta(HiC,len_chr);//generating theta//test
	thEta2 =readTheta(HiC2,len_chr2);//test


	/*for(int i = 0; i < size_b; i++){
                for(int j = 0; j < size_b; j++){
                        string s1 = key_map[i];
                        string s2 = key_map[j];
                        string name = s1 + " " + s2;
                        if(thEta.count(name)){
                                cout<<thEta[s1 + " "+ s2]<<" ";
                        }else{
                                cout<<"0.3 ";
                        }
                }
                cout<<endl;
        }
	exit(1);*/
	//exit(1);
	cout<<"Imputation finished!"<<endl;
	cout<<"Initialization started!"<<endl;
	size_c = 5 * size_b;
	size_c2= 5 *size_b2;
	bool initial_status = false;
	/*while(initial_status == false){
		initial_status = rand_init(size_b, size_c);
	}*/
    std::ifstream file("out8_3");
    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return 1;
    }
	string line;
	while(getline(file, line)){
		vector <double> row;
		istringstream iss(line);
		double value;
		while(iss>>value){
			row.push_back(value);
		}
		slist.push_back(row);
	}
	file.close();
	slist1 = slist;
		
	//mid();
	/*for(int i = 0; i < slist.size(); i++){
		cout<<slist[i][0]<<" "<<slist[i][1]<<" "<<slist[i][1]<<endl;
	}
	exit(1);*/

	//expand_list();
	//exit(1);
	//gen_list(size_b, list_chr_id, list_chr, slist1);
	cout<<"Initialization finished: "<<slist.size()<<" beads"<<endl;	
	//int ck = 38;
	//cout<<cube_io[ck].size()<<" "<<cube_id[ck].size()<<endl;	
	/*for(int i = 0; i < cube_id[ck].size(); i++){
		cout<<cube_id[ck][i]<<" "<<cube_io[ck][i]<<endl;
	}*///check small lists, Verified;
	//start generating small lists for multi-threads
        /*initial_status = false;
        while(initial_status == false){
                initial_status = rand_init(size_b2, size_c2);
        }*/
        //slist = slist2;	
	

	//cout<<"Initialization finished: "<<slist2.size()<<" beads"<<endl;	
	//gen_list(size_b2, list_chr_id2, list_chr2, slist2);
	//exit(1);
	/*for(int i = 0; i < size_b; i++){
		for(int j = 0; j < size_b; j++){
			string s1 = key_map[i];
			string s2 = key_map[j];
			string name = s1 + " " + s2;
			if(HiC.count(name)){
				cout<<HiC[s1 + " "+ s2]<<" ";
			}else{
				cout<<"0 ";
			}
		}
		cout<<endl;
	}*///for ploting the theta matrix;
	
	
	double temperature = 10;
	double tem2 = 10;
	bool cal = false;	
	//cooling(temperature, size_b, map_del, cal);
	//expand_list();
	expand_list();
	cal = true;
	slist = slist2;
	gen_list(size_b2, list_chr_id2, list_chr2, slist2);
	cooling(tem2, size_b2/4, map_del2,cal);	
	/*int counter2 = 0;
	while(counter2<size_c2*10){
		counter2++;
		if(simulate_annealing(10,map_del2,cal)==true){

		}
	}*/
	ofstream fileo;
	fileo.open(out_file);
	/*for(int i = 0; i < slist.size(); i++){
		fileo<<slist[i][0]<<" "<<slist[i][1]<<" "<<slist[i][2]<<endl;
	}*/
	for(int i = 0; i < cube_list.size(); i++){
		for(int j = 0; j < cube_list[i].size();j++){
			fileo<<cube_list[i][j][0]<<" "<<cube_list[i][j][1]<<" "<<cube_list[i][j][2]<<" "<<cube_io[i][j]<<" "<<to_string(int(cube_id[i][j]*res2))<<"-"<<to_string(int((cube_id[i][j]+1)*res2))<<endl;
		}
	}
	
	//fileo.close();
	//cout<<"Modeling finished!!! Input file is "<<hic_file<<" Output coordinates is "<<out_file<<endl;
	/*for(int i = 0; i < size_b; i++){
		cout<<del_map[i]<<endl;
	}*/
        for(int i = 0; i < len_chr2.size(); i++){
                int count = 0;
                for(int j = 0; j < len_chr2[i]; j++){
                        string name = ordered_chr[i] + " " + to_string(j);
                        if(map_del2.count(name)){
                                count++;
                        }
                }
                fileo<<ordered_chr[i]<<" "<<count<<endl;
        }
	fileo.close();
	cout<<"Modeling finished!!! Input file is "<<hic_file<<" Output coordinates is "<<out_file<<endl;
	return 0;
}
