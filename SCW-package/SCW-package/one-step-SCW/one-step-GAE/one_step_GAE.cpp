///remove all zero lines in theta matrix; 
/// and add 0.2tra + 0.8ter, in this version change this, using the count ratio with inter/intra
//direct try 1m for fish
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
int np;
double res;
double res1;
double d0 = 8.0;
int size_b;
double c_ratio;//contacts raiton intra/inter;
vector <int> len_chr;

vector <double> cube_cost;
vector <vector <string>> cube_io;
vector <vector <int>> cube_id;
vector <vector <double>> slist;
vector <vector <vector <double>>> cube_list;

vector <string> ordered_chr;
vector <string> list_chr;
vector <int> list_chr_id;
vector <vector <string>> theta_list;
unordered_map <string, double> thEta;
unordered_map <string, double> HiC;
//unordered_map <string, int> map_key;
unordered_map <string, int> map_del;
unordered_map <int, string> key_map;
unordered_map <int, string> del_map;


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




void readHic(vector<vector <string>> mat){
	int count1 = 0;
	int count2 = 0;
	for(int i = 0; i < mat.size(); i++){
		int s1 = int(stoi(mat[i][1])/res1);
		int s2 = int(stoi(mat[i][3])/res1);
		string name = mat[i][0]+" "+to_string(s1)+" " + mat[i][2] + " " + to_string(s2);
		string name2 = mat[i][2]+" "+to_string(s2)+" " + mat[i][0] + " " + to_string(s1);
		//cout<<name<<endl;
		//cout<<name2<<endl;
		if(mat[i][0] == mat[i][2]){
			count1++;
		}else{
			count2++;
		}
		HiC[name] = 1;
		HiC[name2] = 1;
	}
	cout<<"Inter contacts: "<<count2<<"; Intra contacts: "<<count1<<endl;
	c_ratio = int(count1/count2);
	c_ratio = 1/(c_ratio + 1);
}

void readThe(vector<vector <string>> mat){
	for(int i = 0; i < mat.size(); i++){
                int s1 = int(stoi(mat[i][1]));
                int s2 = int(stoi(mat[i][3]));
                string name = mat[i][0]+" "+to_string(s1)+" " + mat[i][2] + " " + to_string(s2);
                //string name2 = mat[i][2]+" "+to_string(s2)+" " + mat[i][0] + " " + to_string(s1);
                if(stof(mat[i][4]) >= 0.7){
                        thEta[name] = stof(mat[i][4]);
                }
        }

}

void calTheta(vector<vector <int>> ones, int m, int n,int count){
	//vector <vector <int>> ones;
	//cout<<m<<" "<<n<<" "<<ones.size()<<endl;

	for(int i = 0; i < len_chr[m]; i++){
                for(int j = 0; j < len_chr[n]; j++){
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


void readTheta(){
	const size_t THREAD_COUNT = np;//number of thread
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
     
	        	for(int i1 = 0; i1 < len_chr[i]; i1++){
         	        	for(int j1 = 0; j1 < len_chr[j]; j1++){
                	        	string name = ordered_chr[i]+" "+to_string(i1)+" "+ordered_chr[j]+" "+to_string(j1);
                        		//cout<<name<<" "<<HiC[name]<<endl;
					if(HiC.count(name)){
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
			threadPool.emplace_back([ones,i,j,count](){
				calTheta(ones,i,j,count);				
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
			thEta[name] = stod(row[4]);
                }
        }


		
	
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















vector <vector <double>> find_ninit(vector <vector <double>> lis){
	//cout<<lis.size()<<endl;
	vector <vector <double>> possible_pos;
	double x_bottom_bound =((double)lis[lis.size()-1][0])-sqrt(10)-1;
	double x_upper_bound = ((double)lis[lis.size()-1][0])+sqrt(10)+1;
	for(int x = (int)x_bottom_bound; x <= (int)x_upper_bound;x++){
		if(x<0||x>=size_c){
			continue;
		}
		double y_bottom_bound =((double)lis[lis.size()-1][1])-sqrt(10)-1;
        	double y_upper_bound = ((double)lis[lis.size()-1][1])+sqrt(10)+1;
		for(int y = (int)y_bottom_bound; y <= (int)y_upper_bound;y++){
			if(y<0||y>=size_c){
				continue;
			}
			double z_bottom_bound =((double)lis[lis.size()-1][2])-sqrt(10)-1;
                	double z_upper_bound = ((double)lis[lis.size()-1][2])+sqrt(10)+1;
			for(int z = (int)z_bottom_bound; z <= (int)z_upper_bound;z++){
				if(z<0||z>=size_c){
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




bool rand_init(){
        for(int i = 0; i < size_b; i++){
                if(i == 0){
                        timeval t1;
                        gettimeofday(&t1, NULL);
                        srand(t1.tv_usec * t1.tv_sec);

                        int rand1 = rand() % size_c;
                        int rand2 = rand() % size_c;
                        int rand3 = rand() % size_c;
                        vector <double> sl;
                        sl.push_back(rand1);
                        sl.push_back(rand2);
                        sl.push_back(rand3);
                        slist.push_back(sl);
                }else{
                        vector <vector <double>> candidates = find_ninit(slist);

                        if(candidates.size() == 0){
                                cout<<"The size of candidate nodes are 0! initialization failed! There are no enough space. Erase everything and restar"<<endl;
                                slist.clear();
                                return false;
                        }

                        slist.push_back(candidates[0]);
                }
	}
        return true;
}
void mid(){
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

}
void gen_list(){
	
	int len2 = size_b/np + 1;
	int cube_num = int(size_b/len2) + 1;
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
			len3 = size_b - (cube_num-1)*len2;
		}
		vector <vector <double>> clist;
		vector <int> cid;
		vector <string> cio; 
		for(int j = 0; j < len3; j++){
			int k = i*len2 + j;
			//cout<<i<<" "<<k<<"; Chr"<<list_chr[k]<<" : "<<list_chr_id[k]<<endl;
			cid.push_back(list_chr_id[k]);
			cio.push_back(list_chr[k]);
			clist.push_back(slist[k]);
		}
		cube_id[i] = cid;
		cube_io[i] = cio;
		cube_list[i] = clist;
		cube_cost.push_back(0);
	}
	//cout<<len2<<endl;
}


void cal_cube(int num, int m, int n){
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
		if(thEta.count(name)){
			if(thEta[name] > 0.95){
				//if(dist < d0){	
					cost1_1 = cost1_1 + (((dist - d0) * (dist - d0))/(d0*d0));//this 2 is target values,
				//}
			}
			if(thEta[name] >= 0.7 && thEta[name] < 0.95){
				double t1 = thEta[name];
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
	

void call_thread(int m, int n){
	const size_t THREAD_COUNT = cube_cost.size();
	vector <thread> threadPool;
	for(int i = 0; i < cube_cost.size(); i++){
		threadPool.emplace_back([i,m,n](){
			cal_cube(i,m,n);
		});
	}
	for(auto& t: threadPool){
                t.join();
        }
}


double calculate_loss(int m, int n){
	double cost = 0.0;
	for(int i = 0; i < cube_cost.size(); i++){
		cube_cost[i] = 0;
	}
	call_thread(m, n);
	for(int i = 0; i < cube_cost.size(); i++){
		cost = cost + cube_cost[i];
	}
	return cost;
}


bool simulate_annealing(double temperature){
	timeval t1;
	gettimeofday(&t1, NULL);
	srand(t1.tv_usec * t1.tv_sec);
	int rand1 = rand()%cube_list.size();
	int rand2 = rand()%cube_list[rand1].size();
	double cost1 = calculate_loss(rand1, rand2);
	//cout<<cost1<<endl;	
	double oldx = cube_list[rand1][rand2][0];
	double oldy = cube_list[rand1][rand2][1];
	double oldz = cube_list[rand1][rand2][2];
	//cout<<"Old: "<<oldx<<" "<<oldy<<" "<<oldz<<" Cost1: "<<cost1<<endl;
	string name = cube_io[rand1][rand2] + " " + to_string(cube_id[rand1][rand2]);
	int cur = map_del[name];
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

	double cost2 = calculate_loss(rand1,rand2);
	//cout<<"rand cands has "<<cand1.size()<<" Picked "<<rand3<<" th."<<endl;
	//cout<<"New: "<<newx<<" "<<newy<<" "<<newz<<": cost2 is "<<cost2<<endl;

	double Z = cost2 - cost1;
	if(Z <= 0){
		//cout<<" Z <= 0 ture"<<endl<<endl;
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
			return true;
		}else{
			cube_list[rand1][rand2][0] = oldx;
		        cube_list[rand1][rand2][1] = oldy;
		        cube_list[rand1][rand2][2] = oldz;
			slist[cur][0] = oldx;
			slist[cur][1] = oldy;
			slist[cur][2] = oldz;

                        //cout<<"random number is "<<rand4<<" temperate is "<<temperature;
                        //cout<<" check with "<<exp(-Z/temperature)<<" ; large false"<<endl<<endl;
			return false;
		}

	}
}





void cooling(){
	cout<<"Cooling started."<<endl;
	double temperature = 10;
	int failed_times = 0;
	for(int i = 1; i < 1000000; i++){
		cout<<"Temperature: "<<temperature<<endl;
		int counter = 0;
		int accept_counter = 0;
		while(counter < (size_b * 100)){
			counter++;
			if(simulate_annealing(temperature) == true){
				accept_counter++;
			}
			if(accept_counter == (size_b * 10)){
				failed_times = 0;
				break;
			}
			if(counter == ((size_b * 100)-1)){
				failed_times++;
			}
		}
		//cout<<"Accept num: "<<accept_counter<<endl;
		temperature = temperature * 0.9;
		if(failed_times == 3){
			break;
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////int main //////////////////////
int main(int argc, char* argv[]){
	char* hic_file;
	char* length_file;
	char* out_file;
	char* imp_file;
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
		}else if(strcmp(argv[i],"-imp")==0){
			imp_file = argv[i+1];		
		}
	}
	

	res1 = res * 1000000;


	vector <vector <string>> clines = readMatrixFile(length_file);
	cout<<"Chromosome numbers: "<<clines.size()<<endl;
	//int whole_beads = 0;
	int count_beads = 0;// hash{count} = X chr_index
	for(int i = 0; i < clines.size(); i++){
		//cout<<clines[i][1]<<endl;

		//for(int j = 0; j < clines[i].size(); j++){
		//cout<<clines[i][0]<<" "<<clines[i][1]<<endl;
		int num = stoi(clines[i][1]);
		if(num % int(res1) == 0){
			num = num/(res1);
		}else{
			num  = num/(res1) + 1;
		}
		for(int j = 0; j < num; j++){
			string name = clines[i][0] + " " + to_string(j);
			//list_chr.push_back(clines[i][0]);       this will be replaced as bead index: chr;
			//list_chr_id.push_back(j);after editing, this will be replaced as bead index : pos;
			//map_key[name] = count_beads;after this delete
			key_map[count_beads] = name;
			//cout<<name<<" "<<count_beads<<endl;
			count_beads++;
		}
		size_b += num;
		//cout<<clines[i][0]<<" "<<num<<endl;
		len_chr.push_back(num);
		ordered_chr.push_back(clines[i][0]);
		
		//}
	}
	cout<<"Beads numbers: "<<size_b<<endl;
	cout<<"Resolution: "<<res<<" Mb"<<endl;
	for(int i = 0; i < clines.size(); i++){
		cout<<"Chr"<<ordered_chr[i]<<" #beads: "<<len_chr[i]<<endl;
		//cout<<ordered_chr[i]<<" "<<len_chr[i]<<endl; 
	}
	//exit(1);
	
	vector <vector <string> > hicmat = readMatrixFile(hic_file); 
	cout<<"Contacts: "<<hicmat.size()<<endl;
	cout<<"Imputation started!"<<endl;
	readHic(hicmat);// if met an error of segment core, that is the 0 base exceed problem.
	int lt = 0;
	for(int i = 0; i < size_b; i++){
        	int ct = 0; 
		for(int j = 0; j < size_b; j++){
                        string s1 = key_map[i];
                        string s2 = key_map[j];
                        string name = s1 + " " + s2;
                        if(HiC.count(name)){
                                //cout<<HiC[s1 + " "+ s2]<<" ";
                        }else{
                                //cout<<"0 ";
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
			//cout<<"Beads "<<lt<<" is Chr"<<row[0]<<": "<<stoi(row[1])<<endl;
			lt++; 
		}
		if(ct == size_b){
			//cout<<"bead "<<key_map[i]<<"; line "<<i<<" is that line"<<endl;
		}
                //cout<<endl;
        }

	//exit(1);
	size_b = lt;
	vector <vector <string> > impmat = readMatrixFile(imp_file);
	readThe(impmat);//generating theta
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
	cout<<"Imputation finished!"<<endl;
	cout<<"Initialization started!"<<endl;
	size_c = 5 * size_b;
	bool initial_status = false;
	while(initial_status == false){
		initial_status = rand_init();
	}
	mid();
	/*for(int i = 0; i < slist.size(); i++){
		cout<<slist[i][0]<<" "<<slist[i][1]<<" "<<slist[i][1]<<endl;
	}
	exit(1);*/
	gen_list();
	
	//int ck = 38;
	//cout<<cube_io[ck].size()<<" "<<cube_id[ck].size()<<endl;	
	/*for(int i = 0; i < cube_id[ck].size(); i++){
		cout<<cube_id[ck][i]<<" "<<cube_io[ck][i]<<endl;
	}*///check small lists, Verified;
	//start generating small lists for multi-threads
	
	

	cout<<"Initialization finished: "<<slist.size()<<" beads"<<endl;	
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
	
	
	
	cooling();
	
	ofstream fileo;
	fileo.open(out_file);
	for(int i = 0; i < cube_list.size(); i++){
		for(int j = 0; j < cube_list[i].size();j++){
			fileo<<cube_list[i][j][0]<<" "<<cube_list[i][j][1]<<" "<<cube_list[i][j][2]<<" "<<cube_io[i][j]<<" "<<to_string(int(cube_id[i][j]*res1))<<"-"<<to_string(int((cube_id[i][j]+1)*res1))<<endl;
		}
	}
	
	//fileo.close();
	//cout<<"Modeling finished!!! Input file is "<<hic_file<<" Output coordinates is "<<out_file<<endl;
	/*for(int i = 0; i < size_b; i++){
		cout<<del_map[i]<<endl;
	}*/
        for(int i = 0; i < len_chr.size(); i++){
                int count = 0;
                for(int j = 0; j < len_chr[i]; j++){
                        string name = ordered_chr[i] + " " + to_string(j);
                        if(map_del.count(name)){
                                count++;
                        }
                }
                fileo<<ordered_chr[i]<<" "<<count<<endl;
        }
	fileo.close();
	cout<<"Modeling finished!!! Input file is "<<hic_file<<" Output coordinates is "<<out_file<<endl;
	return 0;
}
