#include <sstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h> 
#include <unistd.h>
#include <iomanip>
#include <string.h>
#include <thread>
#include <mutex>
#include <future>
#include <algorithm>
//#include "ThreadPool.h"
//Debug
//#define USEDEBUG
#define PI 3.14159265

#ifdef USEDEBUG
#define Debug(x) std::cout << x
#else
#define Debug(x) 
#endif 

//Debug2
//#define USEDEBUG2

#ifdef USEDEBUG2
#define Debug2(x) std::cout << x
#else
#define Debug2(x)
#endif



using namespace std;
std::mutex mu;
vector <int> XX;
vector <vector <int> > NODES_in1;
vector <vector <double> > BONE;

struct node{
	double x;
        double y;
      	double z;
	double x2 = 0.0;
     	struct node *next;
     	struct node *prev;
}*start;

class Double_list{
	struct node{
        	double x;
        	double y;
        	double z;
		double x2;
        	struct node *next;
        	struct node *prev;
	}*start;

	public:
		int num_nodes;
	public:
		Double_list();
		void append_node(double x, double y, double z);
		void append_x2(double position, double x2);
		double get_x2(double position);
		double get_x_value(double position);
		double get_y_value(double position); //position is starting from 1, end at sequence length, so the first node has 1
		double get_z_value(double position);
		void set_value(double pos, double x, double y, double z);	//position starting from 1
		void clear();	//erase everything, free memory
};

Double_list::Double_list(){
	this -> num_nodes = 0;
	start = NULL;
}

void Double_list::append_node(double x, double y, double z){
	struct node *temp, *current;
	temp = new(struct node);
	temp -> x = x;
	temp -> y = y;
	temp -> z = z;
	temp -> next = NULL;
	if(start == NULL){	//empty list
		temp -> prev = NULL;
		start = temp;
	}
	else{
		current = start;
		while(current -> next != NULL){
			current = current -> next;
		}
		current -> next = temp;
		temp -> prev = current;
	}
	this -> num_nodes++;
}

double Double_list::get_x2(double pos){
        if(start == NULL){      //that means it is empty
                cout << "Fatal Error: Cannot get the value of node, list is empty!\n";
                exit(1);
        }
        if(pos <= 0){
                cout << "Position of the list starting from 1. Invalid position number!\n";
                exit(1);
        }
        struct node *current;
        current = start;
        int i = 1;
        while(i < pos){
                current = current -> next;
                if(current == NULL){
                        cout << "This is beyond the number of node in the list! Cannot get value!\n";
                        exit(1);
                }
                i++;
        }
        return current -> x2;
}


double Double_list::get_x_value(double pos){
	if(start == NULL){	//that means it is empty
		cout << "Fatal Error: Cannot get the value of node, list is empty!\n";
		exit(1);
	}
	if(pos <= 0){
		cout << "Position of the list starting from 1. Invalid position number!\n";
		exit(1);
	}
	struct node *current;
	current = start;
	int i = 1;
	while(i < pos){
		current = current -> next;
                if(current == NULL){
                        cout << "This is beyond the number of node in the list! Cannot get value!\n";
                        exit(1);
                }
		i++;
	}
	return current -> x;
}

double Double_list::get_y_value(double pos){
        if(start == NULL){      //that means it is empty
                cout << "Cannot get the value of node, list is empty!\n";
                exit(1);
        }
        if(pos <= 0){
                cout << "Position of the list starting from 1. Invalid position number!\n";
	          exit(1);
        }
        struct node *current;
        current = start;
        int i = 1;
        while(i < pos){
                current = current -> next;
                if(current == NULL){
                        cout << "This is beyond the number of node in the list! Cannot get value!\n";
                        exit(1);
                }
		i++;
        }
        return current -> y;
}

double Double_list::get_z_value(double pos){
        if(start == NULL){      //that means it is empty
                cout << "Cannot get the value of node, list is empty!\n";
                exit(1);
        }
        if(pos <= 0){
                cout << "Position of the list starting from 1. Invalid position number!\n";
                exit(1);
        }
        struct node *current;
        current = start;
        int i = 1;
        while(i < pos){
                current = current -> next;
		if(current == NULL){
				cout << "This is beyond the number of node in the list! Cannot get value!\n";
				exit(1);
			}
			i++;
		}
		return current -> z;
	}

	void Double_list::set_value(double pos, double x, double y, double z){
		if(start == NULL){      //that means it is empty
			cout << "Cannot set the value of node, list is empty!\n";
			exit(1);
		}
		else{
			struct node *current;
			current = start;
			int i = 1;
			while(i < pos){
				current = current -> next;
				if(current == NULL){
					cout << "This is beyond the number of node in the list! Cannot set value!\n";
					exit(1);
				}
				i++;
			}
			current -> x = x;
			current -> y = y;
			current -> z = z;
		}
	}
	void Double_list::append_x2(double pos, double x2){
                if(start == NULL){      //that means it is empty
                        cout << "Cannot set the value of node, list is empty!\n";
                        exit(1);
                }
		else{
			struct node *current;
			current = start;
			int i = 1;
			while(i < pos){
				current = current -> next;
                                if(current == NULL){
                                        cout << "This is beyond the number of node in the list! Cannot set value!\n";
                                        exit(1);
                                }
                                i++;
			}
			current -> x2 = x2;
		}		
	}
	void Double_list::clear(){
		struct node *pDel = start;
		while(pDel != NULL){
			start = start -> next;
			delete pDel;
			pDel = start;
		}
		start = NULL;
	}


	//The following is a class that performs the simulation
	//
	//
	class Simulation{
			int length_seq2;
			double delta0;
			double theta1;
			double mu1;
			double mu2;
			double rho;
			double beta;	
			double phi;
			double d0;
			double tau;
			double e1;
			double d1;
			double min_dist, max_dist;	//The min and max distance allowed, usually they are [2, sqrt(10)];
		public:
			vector < vector < vector <double> > > list_chr;
			int size_cube;
			int length_seq;
			double temperature;
			int np;
			vector <vector <double> > Ls;
			vector <vector <vector <double> > > hic_list;
			vector <vector <vector <double> > > delta_list;
			vector <vector <vector <int> > > chr_cubes; 
			Double_list list;
			Double_list list2;	
			Simulation(vector <vector <double> > nchr_p, int size, int length, double min, double max, vector < vector < vector <double> > > hic_list,double temper, vector < vector < vector <double> > > delta_list, double r, double the_1, double m_1, double m_2, double tau_1,double b,double del_0,double d_0,double p_1,/*double open_c,*/double d_1,int length2, vector < vector <double> > real_dist2, vector < vector <double> > delta2,double e1,int np);
			bool rand_initialize();	//initialize the sequence in the cube, randomly
			vector<node> find_neighbours_init(Double_list list);	//Find all the possible positions legal. This is for the initialization part,
										//which means the legal positions should be [2, sqrt root of 10] away from all the previous nodes
										//that have existent in the list, which were saved in the list
			void print_pdb_format();
			double calculate_loss(double rand1, double rand2);
			double calculate_cost1(double rand1,double zero,vector< vector <double> > delta,vector < vector <double> > real_dist2,double length_seq);
			double calculate_cost2(int index, int rand);// index is the chr number, rand is the nodes number of the chr
			bool simulate_once(int lstart,int open);	//This function randomly select one node, move it, and make decision of acceptance or reject
			bool simulate_lite(int index, double temperature1);
			void decrease(int index);//doing simulated_annealing jobs
			void multiple();
			vector <vector < vector <int> > > cubes;
			//void test()
			bool simulate_n(int rand0, int cube_unit);
			void print_list_to_file(char *file_name);	//Print final structure into file
			double En = 0.0;
			vector <vector <double> > filter(int r_length);		
			bool expand_list(double res);//expand_list fuction
			double expand_length;
			int cocount = 0;
			bool simulate_cube(double temperature1);
			void call_thread(vector <int> x, vector <vector <int> > nodes_in1, vector <vector <double> > b1, int opz);
			void call_thread2(int chr_size);
			void call_thread4(double rand1, double rand2);
			void call_thread3(int rand1);
			void chr_sum(int rand1, int index, char *file_out);
			void chr_list_sum(int index, int cube_n, int rand1, char *file_name);
			void cal_cube(int cube_num, int rand1, int rand2, int zero);
			vector <vector <vector <int> > > hic_cube;
			void thread_p();
			vector <vector <double> > list_chr_cube_cost;
			vector <double> cube_cost;
			vector <double> rearrenge_cost;
			vector <vector <int> > cnull_line;
			void simulated_annealing(double temperature);
			vector <vector <double> > delta2;
			vector <vector <double> > real_dist;
			vector <vector <double> > real_dist2;
			vector <vector <double> > nchr_p;
			int failed;
			int biggest_length;
			void move_boundary(int r_length);
			void no_move(int r_length, double temperature);
			vector < vector <node> > online_list;// update this list of coordinates
                	vector < vector <string> > olist;
			vector < vector <int> > plist;
			vector < vector <int> > list_start_end;
			double calculate_loss1(double rand1);
			//int index1 = 0;
			//vector <int> nodes_in2;
			//vector <double> b2;
			vector <node> arounding(int index, int random);
			vector<node> around_pos(int zero,int rand1,int length_seq,double d_1);//to find the candidates surrounding the randomly number
			/*std::thread spawn(int index, vector <int> nodes, vector <double> boundary, int zero){
				//mu.lock();
				//return std::thread([=] {simulated_annealing(index, nodes, boundary, zero);});
				//mu.unlock();
			}*/

	};

	Simulation::Simulation(vector <vector <double> > nchr_p, int size, int length, double min, double max, vector < vector < vector <double> > >hic_list, double temper, vector < vector < vector <double> > > delta_list, double r, double the_1, double m_1, double m_2, double tau_1,double b,double del_0,double d_0,double p_1,double d_1,int length2,vector < vector <double> > real_dist2, vector < vector <double> > delta2,double e1,int np){
	this -> hic_cube = hic_cube;
	this -> list_start_end = list_start_end;
	this -> rearrenge_cost = rearrenge_cost;
	this -> list_chr_cube_cost = list_chr_cube_cost;
	this -> cube_cost = cube_cost;
	this -> chr_cubes = chr_cubes;
	this -> size_cube = size;
	this -> length_seq = length;
	this -> length_seq2 = length2;
	this -> min_dist = min;
	this -> max_dist = max;
	this -> hic_list = hic_list;
	this -> real_dist = real_dist;
	this -> real_dist2 = real_dist2;
	this -> delta2 = delta2; 
	this -> delta_list = delta_list;
	this -> temperature = temper;
	this -> rho = r;
	this -> beta = b;
	this -> delta0 = del_0;
	this -> d0 = d_0;
	this -> mu1 = m_1;
	this -> mu2 = m_2;
	this -> tau = tau_1;
	this -> phi = p_1;
	this -> theta1 = the_1;
	this -> d1 = d_1;
	this -> e1 = e1;
	this -> np = np;
	this -> cubes = cubes;
	this -> failed = failed;
}

//For the first node, randomly put it in the cube. For the others call function find_neigbbours_init
//Return true if initilization is successful, i.e., all the nodes have been put in tube, with
//enough space

bool Simulation::rand_initialize(){
	//debug3 cout << "Create a new list of nodes\n";
	Double_list list;
	this -> list = list;
	//cout << "Initialization: randomly put sequence in the cube:\n\n";
	for(int i = 0; i < length_seq; i++){
		//cout<<"	Node number: " << i + 1 << "\n";//@@@
		if(i == 0){	//This is the first node in the list or sequence
			timeval t1;
		        gettimeofday(&t1, NULL);
        		srand(t1.tv_usec * t1.tv_sec);
			//srand (time(NULL));
			int rand1 = rand() % size_cube; //The range of this random number is from 0 to size_cube - 1
                        int rand2 = rand() % size_cube;
                        int rand3 = rand() % size_cube;
			//cout<<"		node put at coordinates: " << rand1 << " " << rand2 << " " << rand3 << "\n";//@@@
			this -> list.append_node(rand1, rand2, rand3);
			////this -> list.append_node(0, 0, 0);	//This doesn't help. trying to maximize the radius of gyration(dont need this, bug sloved, BUG: the new node doesn't have
								//to be smaller than sqrt(10), but just > 2 with each of the nodes already in the list)
								//In this way, every random model will be starting from the origin. This may solve the problem 
								//that there are no enought space as a random first node may leave space not enough to have 
								//available candidate positions satisfying [2, sqrt(10)]
			double n = this -> list.get_x_value(1);
			//cout<<this -> list.get_x_value(1) << this -> list.get_y_value(1) << this -> list.get_z_value(1);//@@@
			//p
			//sleep(1);//@@@	
		}
		else{
			vector<node> candidates = find_neighbours_init(this -> list);
			if(candidates.size() == 0){
				Debug("The size of candidate nodes are 0! initialization failed!There are no enough space. Erase everything and restart\n");
				cout<<"The size of candidate nodes are 0! initialization failed! There are no enough space. Erase everything and restar"<<endl;
				this -> list.clear();
				//sleep(5);
				return false;
			}
			this -> list.append_node(candidates[0].x, candidates[0].y, candidates[0].z);		
			//cout<<candidates[0].x << " " << candidates[0].y << " " << candidates[0].z << " added to list\n";//@@@
			
		}
	}
	 //cout << "Initialization done! Initial nodes are:\n";
	//print_pdb_format();
	return true;
}


void Simulation::print_pdb_format(){
        for(int i = 0; i < length_seq; i++){
		//cout << "ATOM         N   MET A 108                                                   \n";
                 cout << "ATOM  " << right << setw(5) << "    " << " " << setw(4) << "CA" << " " << setw(3) << "MET" << " " << "A" << "        " << setw(8) <<this -> list.get_x_value(i + 1) << setw(8) << this -> list.get_y_value(i + 1) << setw(8) << this ->list.get_z_value(i + 1) << "\n";
		//cout << "ATOM         C   MET A 108                                                   \n";
		//cout << "ATOM         O   MET A 108                                                   \n";
		//cout << "ATOM         CB  MET A 108                                                   \n"; 
		//cout << "ATOM         CG  MET A 108                                                   \n";  
		//cout << "ATOM         SD  MET A 108                                                   \n";
		//cout << "ATOM         CE  MET A 108                                                   \n";
		//cout << "ATOM         H   MET A 108                                                   \n";
        }	
}



void Simulation::print_list_to_file(char *file_name){//use to print the first list
	ofstream file;
	file.open(file_name);
	if(file.is_open()){
		//cout<<"length of the list is" << length_seq << "\n";//@@@
		for(int i = 0; i < length_seq; i++){
			file << setw(8) << this -> list.get_x_value(i + 1) << setw(8) << this -> list.get_y_value(i + 1) << setw(8) << this -> list.get_z_value(i + 1) << "\n"; 
		}
		file.close();
	}
	else{
		cout
 << "Unable to open the output file! Print the coordinates onto screen!";
                for(int i = 0; i < length_seq; i++){
                        cout << setw(8) << this -> list.get_x_value(i + 1) << setw(8) << this -> list.get_y_value(i + 1) << setw(8) << this -> list.get_z_value(i + 1);
                }

	}
}

//Input is a linked list containing all the previously existent nodes. This function at first
//find all the nodes around the last node in the list, i.e., [x-sqrt(10), y - sqrt(10), z - sqrt(10)], x, y, z
//are the coordinates for the last node in the list, and iterate each one of them to 
//see whether it is [2, sqrt(10)] away from the existent nodes. At last, it returns all the legal candidate positions 
vector<node> Simulation::find_neighbours_init(Double_list list){
	vector<node> possible_pos;
	if(list.num_nodes == 0){
		cout << "The doubly linked list is empty, cannot find legal positions for next node!";
		exit(1);
	}
	//This newly added nodes must be within [2, sqrt(10)] of the last node in the list
	//so the new node must locate in [x - max - 1, x + max + 1], +-1 is because rounding, for
	//safe reason, because C truncate double when convert to int, so include more possibilities
	//So get all of the possible positions based on this criteria first, then graduately
	//eliminate unqualified ones
	double x_bottom_bound = ( (double)list.get_x_value(list.num_nodes) ) - (this -> max_dist) - 1;
	double x_upper_bound = ( (double)list.get_x_value(list.num_nodes) ) + this -> max_dist + 1;
	for(int x = (int)x_bottom_bound; x <= (int)x_upper_bound; x++){
		if(x < 0 || x >= this->size_cube){	//beyond the boundery of the cube, cube coordinates starting from 0
			continue;
		}
//		cout<<"x is " << x << "\n";//@@@
		double y_bottom_bound = ( (double)list.get_y_value(list.num_nodes) ) - this -> max_dist - 1;
		double y_upper_bound = ( (double)list.get_y_value(list.num_nodes) ) + this -> max_dist + 1;
		for(int y = (int)y_bottom_bound; y <= (int)y_upper_bound; y++){
			if(y < 0 || y >= this -> size_cube){
				continue;
			}
//			cout<<"y is " << y << "\n";//@@@
			double z_bottom_bound = ( (double)list.get_z_value(list.num_nodes) ) - this -> max_dist - 1;
			double z_upper_bound = ( (double)list.get_z_value(list.num_nodes) ) + this -> max_dist + 1;
			for(int z = (int)z_bottom_bound; z <= (int)z_upper_bound; z++){
				if(z < 0 || z >= this -> size_cube){
					continue;
				}
//				cout<<"z is " << z << "\n";//@@@
				struct node *temp;
			        temp = new(struct node);
			        temp -> x = x;
			        temp -> y = y;
			        temp -> z = z;
				//cout<<"		" << x << " " << y << " " << z << " " << "added as candidate positions\n";//@@@
				possible_pos.push_back(*temp);
				delete temp;
			}
		}
	}
//	cout<<"		number of candidates now: " << possible_pos.size() << "\n";//@@@

	//These possible positions should not contain any ones that are without [2, sqrt(10)] towards
	//to any of the pre-existed nodes saved in list. Now get ride of these ndoes.
	//cout<<"		Checking distance between candidate positions and each node in list\n";//@@@
	vector<node> possible_pos_2;	//Define another vector saving the legal positions
	bool random_find=false;
	int count=0;
	while(random_find==false){
		count++;
		if(count==possible_pos.size()){
			break;
		}
		timeval t1;
                gettimeofday(&t1, NULL);
                srand(t1.tv_usec * t1.tv_sec);
		
		int rand1= rand() % possible_pos.size()-1;	
		int x_p = possible_pos[rand1].x;
		int y_p = possible_pos[rand1].y;
		int z_p = possible_pos[rand1].z;
		//cout<<"Checking candidate position with each nodes already in the list..." << m << "\n";//@@@
		//cout<<x_p << " " << y_p << " " << z_p << "\n"; //@@@
		bool removed = false;
		for(int i = 1; i <= list.num_nodes; i++){	//List position starting from 1
			//cout<<"		checking the " << i << " node in the list already\n";	//@@@
			double x_e = list.get_x_value(i);
			double y_e = list.get_y_value(i);
			double z_e = list.get_z_value(i);
			//cout<<"		calculating the dist to node " << x_e << " " << y_e << " " << z_e << "\n";//@@@
			double dist = sqrt( (x_p - x_e) * (x_p - x_e) + (y_p - y_e) * (y_p - y_e) + (z_p - z_e) * (z_p - z_e) );
			//cout<<"		distance = " << dist << "\n";//@@@
			if(i < list.num_nodes){
				if(dist < this -> min_dist || dist == sqrt(8)){		//|| dist > this -> max_dist is removed, because it doesn't need to be smaller than sqrt(10) with all the other nodes, but only 
										//the last node in the list
					removed = true;
					//cout<<"		removed\n";//@@@
					break;
				}
			}
			if(i == list.num_nodes){	//This is judging with the last node in the list, this one can use dist > this -> max_dist
				if(dist < this -> min_dist || dist == sqrt(8) || dist > this -> max_dist){
					removed = true;
					//cout<<"         removed\n";//@@@
					break;
				}
			}
		}
		if(removed == false){
			possible_pos_2.push_back(possible_pos[rand1]);
			random_find=true;
		}
	}
	if(possible_pos_2.size() < 0){
		cout << "Fatal Error: the number of possible positions cannot be a negative number. Something is very wrong!\n";
		exit(1);
	}
	//cout<<"		size of candidates kept: " << possible_pos_2.size() << "\n";//@@@changed debug
	return possible_pos_2;
}


//Read from a file that contains the correlation of Hi-C, return a vector of vector
vector <string> read_chr_nums(char* file){
	ifstream myfile;
	myfile.open(file);
	string line;
	string temp;
	vector <string> s1;
	while (getline(myfile, line)){//save two columns position in array
                istringstream istr(line);
		vector <string> row;
                while (istr >> temp) {
			string tempf = temp.c_str();
                	row.push_back(temp);
                }
		if(row[0] == row[2]){
			//cout<<row[0]<<" "<<endl;
			int oz = 0;
			for(int i = 0; i < s1.size(); i++){
				if(s1[i] == row[0]){
					oz = 1;
				}
			}
			if(oz == 0){
				s1.push_back(row[0]);
			}
		}
                row.clear();
                istr.clear();
                line.clear();
        }
	/*for(int i = 0; i < s1.size(); i++){
		cout<<s1[i]<<endl;
	}*/
	//cout<<s1.size()<<endl;
	return s1;
	//exit(1);
}
vector <vector <vector <double> > > read_w_hic(double res, char* file, vector <string> cnum, vector <vector <string> > cname){//read whole genom hic
	vector < vector < vector <double> > > w1;
	//vector <double> big;
	for(int i = 0; i < cnum.size(); i++){
	        ifstream myfile;
        	myfile.open(file);
        	string line;
        	string temp;
		double big = 0.0;
		vector <vector <double> > mat;
		while (getline(myfile, line)){
			vector <string> row;
			istringstream istr(line);
			while (istr >> temp){
                        	/*if(temp == "X"){
                                	temp = "20";
                        	}
				if(temp == "Y"){
					temp = "21";
				}*/
				string tempf = temp.c_str();
				/*if(big <= tempf){
					big = tempf;
				}*/
                        	row.push_back(tempf);			
			}
			if(row[0] == row[2] && row[0] == cnum[i]){
				vector <double> row2;
				double temp1 = atof(row[1].c_str());
				double temp3 = atof(row[3].c_str());
				temp1 = temp1 / (res * 1000000) + 1;
				temp3 = temp3 / (res * 1000000) + 1;
				/*if(big <= temp1){
                                	big = temp1;
                                }
				if(big <= temp3){
					big = temp3;
				}*/
				row2.push_back(temp1);
				row2.push_back(temp3);
				mat.push_back(row2);
		
				//cout<<row[0]<<" "<<temp1<<" "<<row[2]<<" "<<temp3<<endl;
			}
		
		}
                int length_c;
                for(int j = 0; j < cname.size(); j++){
                        if(cnum[i] == cname[j][0]){
                                int lcc = atoi(cname[j][1].c_str());
                                double res2 = res * 1000000;
                                length_c = int(lcc/res2) + 1;
                        }
                }

		//exit(1);
		vector <vector <double> > w2;
		for(int j = 0; j < length_c; j++){
			vector <double> row;
			for(int k = 0; k < length_c; k++){
				double temp_n = 0.0;
				row.push_back(temp_n);
			}
			w2.push_back(row);
		}
		for(int j = 0; j < mat.size(); j++){
			double i0 = mat[j][0];
			double i1 = mat[j][1];
			if(i0 != i1){
				w2[i0 - 1][i1 - 1] = 1;
				w2[i1 - 1][i0 - 1] = 1;
			}
		}
		cout<<"chr "<<cnum[i]<<" "<<length_c<<endl;
		w1.push_back(w2);
		//cout<<i<<endl;
				
	}
	return w1;
}
vector < vector <double> > index_l;
/////////////////////
vector< vector<double> > read_correlation(double res, int size, char* file, vector < vector <string> > index){//read two position bead pairs to generate a hic matrix
        vector< vector<string> > vector2;
        vector< vector<double> > vector3;             
	vector<string> row;
        string line;
        string temp;
        ifstream myfile;
        myfile.open(file);
        if (! myfile.is_open()) {
                cout<<"File :< " << file << " >open error"<< endl;
                exit(EXIT_FAILURE);
        }
        int big = 0;
	res = res * 1000000;
        while (getline(myfile, line)){//save two columns position in array
                istringstream istr(line);
                while (istr >> temp) {
			/*if(temp == "X"){
				temp = "20";
			}
			if(temp == "Y"){
				temp = "21";
			}
			double tempf = atof(temp.c_str());*/
			
                        row.push_back(temp);
                        //if(big<=tempf){
                         //       big=tempf;//find the biggest position int the hic contact, set it as the nodes number of the chromosome
                        //}
                }
                //cout<<endl;
                double ar1 = atof(row[1].c_str());
		double ar2 = atof(row[3].c_str());
                ar1 = (ar1 / res) + 1;
		ar2 = (ar2 / res) + 1;
		row[1] = to_string(ar1);
		row[3] = to_string(ar2);
		//cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<" "<<row[3]<<endl;
                vector2.push_back(row);
                if(vector2.back().empty()){
                        cout<<"File :< " << file << " >has an empty row"<< endl;
                        exit(EXIT_FAILURE);
                }
                row.clear();
                istr.clear();
                line.clear();
        }
	//cout<<"!#"<<endl;
        //cout<<"biggest is "<<big<<endl;
        myfile.close();
        //res = res*1000000;
        /*int big_num = big/res+1;
        //cout<<"big bead is "<<big_num<<endl;
        vector <int> node_big;
        for(int i = 0; i < index.size(); i++){
                int k  = 0;
                node_big.push_back(k);
        }*/
	int all_node_beads = 0;
	for(int i = 0; i < index.size(); i++){
		all_node_beads = all_node_beads + atoi(index[i][1].c_str());
	}
	vector <vector <int> > string_compair_int;
	int count = 0;
	for(int i = 0; i < index.size(); i++){
		vector <int> small_int;
		int st = count;
		int st1 = count + atoi(index[i][1].c_str());
		count = count + atoi(index[i][1].c_str()); 
		for(int j = st; j < st1; j++){
			small_int.push_back(j);
		}
		string_compair_int.push_back(small_int);
		cout<<"from "<<st<<":"<<st1<<" ;size "<<small_int.size()<<endl;
	}
	cout<<all_node_beads<<endl;
	vector <vector <double> > node_m;
	for(int i = 0; i < all_node_beads; i++){
                vector <double> l;
                for(int j = 0; j < all_node_beads; j++){
                        l.push_back(0);
                }
                node_m.push_back(l);
        }
	/*
	cout<<"!$"<<endl;
	int num = 0;
	for(int i = 0; i < index.size(); i++){
		for(int j = 0; j < index[i][1]; j++){
			vector <double> l;
			num++;
			l.push_back(i+1);
			l.push_back(j+1);
			l.push_back(num);
			//cout<<i + 1<<" "<<j + 1<<" "<<num<<endl;
			index_l.push_back(l);
		}
	}*/
	//exit(1);
	//vector <vector <int> > pairs;
	cout<<" original vector has "<<vector2.size()<<" pairs"<<endl;
	sleep(1);
	for(int i = 0; i < vector2.size(); i++){
		int x1 = 0;
		int x2 = 0;
		for(int j = 0; j < index.size(); j++){
			if(vector2[i][0] == index[j][0]){
				x1 = string_compair_int[j][atoi(vector2[i][1].c_str()) - 1];
			}
			if(vector2[i][2] == index[j][0]){
				x2 = string_compair_int[j][atoi(vector2[i][3].c_str()) - 1];
			}
		}
		//cout<<x1<<" "<<x2<<" pairs is 1"<<endl;
		node_m[x1][x2] = 1;
		node_m[x2][x1] = 1;	
	}
	//exit(1);
	//cout<<"pairs size "<<pairs.size()<<endl;
	//sleep(1);
	//cout<<"!@"<<endl;
	//vector < vector <double> > node_m;
	/*for(int i = 0; i < all_node_beads; i++){
		vector <double> l;
		for(int j = 0; j < all_node_beads; j++){
			l.push_back(0);
		}
		node_m.push_back(l);
	}*/
	
	//cout<<node_m.size()<<endl;
	//cout<<"!%"<<endl;
	//cout<<pairs.size()<<endl;
	//cout<<"res "<<res<<endl;
	//if(res == 500000){
	//for(int k = 0; k < pairs.size(); k++){
	//	node_m[pairs[k][0] - 1][pairs[k][1] - 1] = 1;
	//	node_m[pairs[k][1] - 1][pairs[k][0] - 1] = 1;
		//cout<<k<<" "<<pairs[k][0]<<" "<<pairs[k][1]<<" : "<<pairs[k][1]<<" "<<pairs[k][0]<<endl;
	//}
	/*for(int i = 0; i < node_m.size(); i++){
		for(int j = 0; j < node_m[i].size(); j++){
			cout<<node_m[i][j]<<" ";
		}
		cout<<endl;
	}*/
	//exit(1);
	return node_m;
}

vector < vector <double> > calculate_t2(double length_seq, vector < vector <double> > real_dist,double mu2,double d0){
        vector< vector<double> > t1;
        vector< vector<double> > t2;
        vector< vector<double> > t3;

        for(int i=0; i<length_seq;i++){//get all the real hic contact bead pairs
                for(int j=0; j<length_seq;j++){
                        if(real_dist[i][j]==1){
                                vector <double> row;
                                row.push_back(i);
                                row.push_back(j);
                                t2.push_back(row);
                        }
                }
        }
	for(int i = 0; i < length_seq; i++){
                vector<double> row;
                for(int j = 0; j<length_seq; j++){
                    double temp1=0;
                    if(j>=i){
                        if(real_dist[i][j]==0 && i != j){
                                for(int k =0; k< t2.size(); k++){
                                        if(t2[k][0]==i && t2[k][1]==j){
                                                temp1=1;
                                        }else{
                                                double i1 = t2[k][0];
                                                double j1 = t2[k][1];
                                                double temp2= (i1-i)*(i1-i)+(j1-j)*(j1-j);
                                                temp2 = temp2/mu2;
                                                temp2 = exp(-temp2);
                                                temp1 = temp1 + temp2;
                                        }
                                }
				if(temp1 > 1){
					temp1 = 1;
                                }
                        }else{
                                temp1=0;
                        }
                    }else{
		    }
		    row.push_back(temp1);	
                }
                t3.push_back(row);
        }

        for(int i=0;i<length_seq;i++){//put the half matrix to a full one
                vector <double> row;
                for(int j=0;j<length_seq;j++){
                        double temp1=0;
                        if(i<=j){
                                temp1=t3[i][j];
                        }
                        if(i>j){
                                temp1=t3[j][i];
                        }
                        row.push_back(temp1);
                }
                t1.push_back(row);
        }
        return t1;
}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculate the delta matrix
vector< vector<double> > calculate_t1(double length_seq, vector < vector <double> > real_dist,double mu2,double d0, vector < vector <string> > index){//to generate the delta matrix, d0 is useless here
        vector< vector<double> > t1;
        vector< vector<double> > t2;
        vector< vector<double> > t3;
	//cout<<"length is "<<length_seq<<endl;
	//exit(1);
        for(int i = 0; i < length_seq; i++){
		vector <double> t11;
		for(int j = 0; j < length_seq; j++){
			t11.push_back(0);
		}
		t1.push_back(t11);
	}
	int all_node_beads = 0;
	for(int i = 0; i < index.size(); i++){
		all_node_beads = all_node_beads + atoi(index[i][1].c_str());
	}
	vector <vector <int> > string_compair_int;
        int count = 0;
        for(int i = 0; i < index.size(); i++){
                vector <int> small_int;
                int st = count;
                int st1 = count + atoi(index[i][1].c_str());
                count = count + atoi(index[i][1].c_str());
                for(int j = st; j < st1; j++){
                        small_int.push_back(j);
                }
                string_compair_int.push_back(small_int);
                //cout<<"from "<<st<<":"<<st1<<" ;size "<<small_int.size()<<endl;
        }

	for(int i = 0; i < index.size(); i++){
		for(int j = i; j < index.size(); j++){
			vector <vector <int> > pairs;
			int st1 = string_compair_int[i].size();
			int st2 = string_compair_int[j].size();
			for(int k = string_compair_int[i][0]; k < string_compair_int[i][0] + st1; k++){
				for(int l = string_compair_int[j][0]; l < string_compair_int[j][0] + st2; l++){
					if(real_dist[k][l] == 1){
						vector <int> row;
						row.push_back(k);
						row.push_back(l);
						pairs.push_back(row);
						//cout<<"has pairs "<<k<<" "<<l<<endl;	
					}
				}
			}
			for(int k = string_compair_int[i][0]; k < string_compair_int[i][0] + st1; k++){
                                for(int l = string_compair_int[j][0]; l < string_compair_int[j][0] + st2; l++){
					double temp1 = 0.0;
					if(real_dist[k][l] == 0){
						for(int a = 0; a < pairs.size(); a++){
							int k1 = pairs[a][0];
							int l1 = pairs[a][1];
							double temp2 = (k1 - k) * (k1 - k) + (l1 - l) * (l1 - l);
							temp2 = temp2 / mu2;
							temp2 = exp(-temp2);
							temp1 = temp1 + temp2;
						}
						if(temp1 > 1){
							temp1 = 1;
						}
					}else{
						temp1 = 1;
					}
					t1[k][l] = temp1;
					t1[l][k] = temp1;
				}
			}
		}
	}
	return t1;
}
int cal_big(vector <double> coord){
	double small = 1000000;
        double big = -1000000;
        for(int i = 0; i < coord.size(); i++){
                if( small > coord[i]){
                        small = coord[i];
                }
                if(big < coord[i]){
                        big = coord[i];
                }
        }
	double length = big - small;	
	return length;
}
int calculate_middle_point1(vector <double> coord){
        double small = 1000000;
        double big = -1000000;
        double sum = 0;
        for(int i = 0; i < coord.size(); i++){
                if( small > coord[i]){
                        small = coord[i];
                }
                if(big < coord[i]){
                        big = coord[i];
                }
                sum = sum + coord[i];
        }
        sum = sum / coord.size();
        double middle = (small + big)/2;
        return sum;
}

vector <node> change_position1(vector <vector <node> > ol, vector <node> v1, int length, int biggilen){//20 should be changed into biggest_length
        int a, b, c, a2, b2, c2;
        vector <node> v2;
        vector <double> c_x;
        vector <double> c_y;
        vector <double> c_z;
        for(int i = 0; i < v1.size(); i++){
                c_x.push_back(v1[i].x);
                c_y.push_back(v1[i].y);
                c_z.push_back(v1[i].z);
        }
        double m_x = calculate_middle_point1(c_x);
        double m_y = calculate_middle_point1(c_y);
        double m_z = calculate_middle_point1(c_z);
        int k = 0;
        while(k == 0){
	        timeval t1;
        	gettimeofday(&t1, NULL);
        	srand(t1.tv_usec * t1.tv_sec);
                int rand1 = rand() % (length * biggilen) + 1;
                int rand2 = rand() % (length * biggilen) + 1;
                int rand3 = rand() % (length * biggilen) + 1;
                double t_x = double(rand1) - m_x;
                double t_y = double(rand2) - m_y;
                double t_z = double(rand3) - m_z;
                for(int i = 0; i < v1.size(); i++){
                        double x = v1[i].x + t_x;
                        double y = v1[i].y + t_y;
                        double z = v1[i].z + t_z;
			double x2 = v1[i].x2;
                        struct node*temp;
                        temp = new(struct node);
                        temp -> x = x;
                        temp -> y = y;
                        temp -> z = z;
			temp -> x2 = x2;
                        v2.push_back(*temp);
                        delete temp;
                }
                int count = 0;
                for(int i = 0; i < v2.size(); i++){
                        for(int j = 0; j < ol.size(); j++){
                                for(int l = 0; l < ol[j].size(); l++){
                                        if(v2[i].x == ol[j][l].x && v2[i].y == ol[j][l].y && v2[i].z == ol[j][l].z && ol[j][l].x2 != 1){
                                                count++;
                                        }
                                        if(v2[i].x < 0 || v2[i].x > (length * biggilen) || v2[i].y < 0 || v2[i].y > (length * biggilen) || v2[i].z < 0 || v2[i].z > (length * biggilen)){
                                                count++;
                                        }
                                }
                        }
                }
                if(count == 0){
                        k = 1;
                }else{
                        v2.clear();
                }
        }
        return v2;

}

vector <node> rotate_coord(vector <node> v1, int length){
        int a = 0;
	int b = 0;
	int c = 0;
	int d = 0;
        vector <node> v2;
        vector <double> oldx;
        vector <double> oldy;
        vector <double> oldz;
	vector <double> oldx2;
        for(int i = 0; i< v1.size(); i++){
                double x = v1[i].x - a*length*5;
                double y = v1[i].y - b*length*5;
                double z = v1[i].z - c*length*5;
		double x2 = v1[i].x2;
                oldx.push_back(x);
                oldy.push_back(y);
                oldz.push_back(z);
		oldz.push_back(x2);
                struct node*temp;
                temp = new(struct node);
                temp -> x = x;
                temp -> y = y;
                temp -> z = z;
		temp -> x2 = x2;
                v2.push_back(*temp);
                delete temp;
	}
        double o_x = calculate_middle_point1(oldx);
        double o_y = calculate_middle_point1(oldy);
        double o_z = calculate_middle_point1(oldz);

        timeval t1;
        gettimeofday(&t1, NULL);
        srand(t1.tv_usec * t1.tv_sec);
        int rand1 = rand() % 360;
        int rand2 = rand() % 360;
        int rand3 = rand() % 360;
        vector <double> arrayx;
        vector <double> arrayy;
        vector <double> arrayz;
        vector <node>  a_r;

        for(int i = 0; i<v2.size(); i++){
                double theta_a = rand1*PI/360;
                double theta_b = rand2*PI/360;
                double theta_c = rand3*PI/360;
                double sin_a = sin(theta_a);
		double cos_a = cos(theta_a);
                double sin_b = sin(theta_b);
                double cos_b = cos(theta_b);
                double sin_c = sin(theta_c);
                double cos_c = cos(theta_c);
                double r_x = (cos_a*cos_c-cos_b*sin_a*sin_c)*v2[i].x + (-cos_b*cos_c*sin_a-cos_a*sin_c)*v2[i].y + (sin_a*sin_b)*v2[i].z;
                double r_y = (cos_c*sin_a+cos_a*cos_b*sin_c)*v2[i].x + (cos_a*cos_b*cos_c-sin_a*sin_c)*v2[i].y + (-cos_a*sin_b)*v2[i].z;
                double r_z = (sin_b*sin_c)*v2[i].x + (cos_c*sin_b)*v2[i].y + cos_b*v2[i].z;

                arrayx.push_back(r_x);
                arrayy.push_back(r_y);
                arrayz.push_back(r_z);
                struct node*temp;
                temp = new(struct node);
                temp -> x = r_x;
                temp -> y = r_y;
                temp -> z = r_z;
                a_r.push_back(*temp);
                delete temp;

        }
        double m_x = calculate_middle_point1(arrayx);
        double m_y = calculate_middle_point1(arrayy);
        double m_z = calculate_middle_point1(arrayz);
        double movex = o_x - m_x;
        double movey = o_y - m_y;
        double movez = o_z - m_z;
        vector <node> v3;
        for(int i = 0; i < a_r.size(); i++){
                double x = a_r[i].x + movex + double(a*length*5);
                double y = a_r[i].y + movey + double(b*length*5);
                double z = a_r[i].z + movez + double(c*length*5);
		double x2 = v1[i].x2;
                struct node*temp;
                temp = new(struct node);
                temp -> x = x;
                temp -> y = y;
                temp -> z = z;
		temp -> x2 = x2;
                v3.push_back(*temp);
                delete temp;
        }
        return v3;
}
double Simulation::calculate_loss1(double rand1){
	double cost = 0.0;
	call_thread3(rand1);
	for(int i = 0; i < rearrenge_cost.size(); i++){
		//cout<<"rearrenge_cost "<<i<<": "<<rearrenge_cost[i]<<" ;"<<endl;
		cost = cost + rearrenge_cost[i];
	}
	return cost;
}
//this function is used to calcualte loss_function
double Simulation::calculate_loss(double rand1, double rand2){//calculate loss function for each small cube
	double cost = 0.0;
	/*cout<<"before multi"<<endl;
	for(int i = 0; i < cube_cost.size(); i++){
		cout<<cube_cost[i]<<" is "<<i+1<<endl;
	}*/
	call_thread4(rand1, rand2);
	//exit(1);
	for(int i = 0; i < cube_cost.size(); i++){
		cost = cost + cube_cost[i];
		//cout<<cube_cost[i]<<endl;
	}
	//cout<<"total cost is "<<cost<<endl;
	//exit(1);
	return cost;	
}
double Simulation::calculate_cost1(double rand1,double zero,vector< vector <double> > delta,vector< vector <double> > real_dist,double length_seq){
	double cost1 = 0.0;
	double cost1_1 = 0.0;
	double cost1_2 = 0.0;
	double cost1_3 = 0.0;
	for(int i = 0; i < length_seq; i++){
		if(rand1 != i + 1){//
	//		cout<<"rand1 is "<<rand1<<";delta0 is "<<delta0<<"; theta1 is "<<theta1<<"; mu1 is "<<mu1<<"; beta is "<<beta<<"; tau is  "<<tau<<"; phi is "<<phi<<"; rho is "<<rho<<";"<<endl;
			double x_i,y_i,z_i,x_j,y_j,z_j;
			if(zero==1){
				x_i = list.get_x_value(i + 1);
        	                y_i = list.get_y_value(i + 1);
                	        z_i = list.get_z_value(i + 1);
                        	x_j = list.get_x_value(rand1);
                        	y_j = list.get_y_value(rand1);
                        	z_j = list.get_z_value(rand1);
			}
			if(zero==2){
	                        x_i = list2.get_x_value(i + 1);
        	                y_i = list2.get_y_value(i + 1);
                	        z_i = list2.get_z_value(i + 1);
                        	x_j = list2.get_x_value(rand1);
                        	y_j = list2.get_y_value(rand1);
                        	z_j = list2.get_z_value(rand1);
			}
                        double dist = sqrt( (x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j) + (z_i - z_j) * (z_i - z_j) );
                        if(real_dist[rand1-1][i]==1 || delta[rand1-1][i]==1){
				cost1_1 = cost1_1 + ((dist-delta0)*(dist-delta0))/(delta0*delta0);
			}
			if(real_dist[rand1-1][i]==0 && delta[rand1-1][i]<1 && delta[rand1-1][i]>theta1){
                                double t1 = delta[rand1-1][i];
                                t1 = pow(t1,1.0/3);
                                t1 = delta0/t1;
                                cost1_2=cost1_2 + beta*(1- exp(-(dist-t1)*(dist-t1)/mu1));
                        }
                        if(real_dist[rand1-1][i]==0 && delta[rand1-1][i]<=theta1){
                                double t2 = 0;
                                t2 = pow(theta1,1.0/3);
                                t2 = delta0/t2;
                                cost1_3=cost1_3 + tau*(1- 1/(1+exp(-(dist-(t2-rho)))/phi)) ;
			}
		}
	}
	cost1=cost1_1+cost1_2+cost1_3;
	//cout<<"cost is "<<cost1<<endl;
	return cost1;
}
void Simulation::chr_list_sum(int index, int cube_n, int rand1, char *file_name){//index is chr num, rand1 is the nodes num of chr
	double cost1 = 0.0;
	double cost1_1 = 0.0;
	double cost1_2 = 0.0;
	double cost1_3 = 0.0;
	//ofstream cost_out;
	//cost_out.open(file_name, fstream::app);
	double x_i = list_chr[index][rand1][0];
        double y_i = list_chr[index][rand1][1];
        double z_i = list_chr[index][rand1][2];
	for(int i = 0; i < chr_cubes[index][cube_n].size(); i++){
        	double x_j = list_chr[index][chr_cubes[index][cube_n][i]][0];
                double y_j = list_chr[index][chr_cubes[index][cube_n][i]][1];
                double z_j = list_chr[index][chr_cubes[index][cube_n][i]][2];
                double dist = sqrt( (x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j) + (z_i - z_j) * (z_i - z_j) );
		//cost_out<<"cube is "<<cube_n<<"; node i is "<<chr_cubes[index][cube_n][i]<<"; dist is "<<dist<<" cost is "<<cost1_1<<" "<<cost1_2<<" "<<cost1_3<<endl;
                if(hic_list[index][rand1][chr_cubes[index][cube_n][i]] == 1 || delta_list[index][rand1][chr_cubes[index][cube_n][i]] == 1){
                	cost1_1 = cost1_1 + ((dist-delta0)*(dist-delta0))/(delta0*delta0);
                }
                if(hic_list[index][rand1][chr_cubes[index][cube_n][i]] == 0 && delta_list[index][rand1][chr_cubes[index][cube_n][i]] < 1 && delta_list[index][rand1][chr_cubes[index][cube_n][i]] > theta1){
                	double t1 = delta_list[index][rand1][chr_cubes[index][cube_n][i]];
                        t1 = pow(t1,1.0/3);
                        t1 = delta0/t1;
                        cost1_2=cost1_2 + beta*(1- exp(-(dist-t1)*(dist-t1)/mu1));
                }
                if(hic_list[index][rand1][chr_cubes[index][cube_n][i]] == 0 && delta_list[index][rand1][chr_cubes[index][cube_n][i]] <= theta1){
                	double t2 = 0;
                        t2 = pow(theta1,1.0/3);
                        t2 = delta0/t2;
                        cost1_3=cost1_3 + tau*(1- 1/(1+exp(-(dist-(t2-rho)))/phi)) ;
                }
		//cost_out<<"cube "<<cube_n<<"; node i is "<<chr_cubes[index][cube_n][i]<<"; cost is "<<cost1_1<<" "<<cost1_2<<" "<<cost1_3<<endl;
	}
	cost1 = cost1_1 + cost1_2 + cost1_3;
	//cost_out.close();
	list_chr_cube_cost[index][cube_n] = cost1;
}
double Simulation::calculate_cost2(int index, int rand1){
        double cost1 = 0.0;
	/*call_thread2(index, rand1);//index of cube;
	for(int i = 0; i < list_chr_cube_cost[index].size(); i++){
		cost1 = cost1 + list_chr_cube_cost[index][i];//should generate a new vector 
	}
	return cost1;*/
        double cost1_1 = 0.0;
        double cost1_2 = 0.0;
        double cost1_3 = 0.0;
        double x_i = list_chr[index][rand1][0];
        double y_i = list_chr[index][rand1][1];
        double z_i = list_chr[index][rand1][2];
	int lenq = list_chr[index].size();
        for(int i = 0; i < lenq; i++){
		if(i != rand1){
                       // cout<<"rand1 is "<<rand1<<";delta0 is "<<delta0<<"; theta1 is "<<theta1<<"; mu1 is "<<mu1<<"; beta is "<<beta<<"; tau is  "<<tau<<"; phi is "<<phi<<"; rho is "<<rho<<";"<<endl;
                        double x_j = list_chr[index][i][0];
                        double y_j = list_chr[index][i][1];
                        double z_j = list_chr[index][i][2];
                        double dist = sqrt( (x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j) + (z_i - z_j) * (z_i - z_j) );
                        if(delta_list[index][rand1][i] >= 0.95 && delta_list[index][rand1][i] < 1){
                                cost1_1 = cost1_1 + ((dist-delta0)*(dist-delta0))/(delta0*delta0);
                        }
                        if(delta_list[index][rand1][i] >= 0.7 && delta_list[index][rand1][i] < 0.95){
                                double t1 = delta_list[index][rand1][i];
                                t1 = pow(t1,1.0/3);
                                t1 = delta0/t1;
                                cost1_2=cost1_2 + beta*(1- exp(-(dist-t1)*(dist-t1)/mu1));
                        }
                        if(delta_list[index][rand1][i] < 0.7){
                                double t2 = 0;
                                t2 = pow(0.7,1.0/3);
                                t2 = delta0/t2;
                                cost1_3=cost1_3 + tau*(1- 1/(1+exp(-(dist-(t2-rho)))/phi)) ;
                        }
                }
        }
        cost1=cost1_1+cost1_2+cost1_3;
	//cout<<"cost is "<<cost1<<endl;
        return cost1;
}
vector <node> Simulation::arounding(int index, int random){
	vector <node> candidates;
	int lenq = list_chr[index].size();
	double x, y , z, xmin, xmax, ymin, ymax, zmin, zmax;
	x = list_chr[index][random][0];
	y = list_chr[index][random][1];
	z = list_chr[index][random][2];
	//cout<<"xyz are "<<x<<" "<<y<<" "<<z<<", random number is "<<random<<"; d1 is "<<d1<<endl;
	double bo = list_chr[index].size();
	for(double i = x - 1; i <= x + 1; i++ ){
		for(double j = y - 1; j <= y + 1; j++){
			for(double k = z - 1; k <= z + 1; k++){
				if(x == i && y == j && z == k){ continue;}
				else{
					//cout<<"ijk are"<<i<<" "<<j<<" "<<k<<endl;
					if(random == 0){
						xmax = list_chr[index][1][0];
						ymax = list_chr[index][1][1];
						zmax = list_chr[index][1][2];
						double dist2 = sqrt((i-xmax)*(i-xmax)+(j-ymax)*(j-ymax)+(k-zmax)*(k-zmax));
						if(dist2 > d1){
							continue;
						}
					}else if(random == lenq - 1){
						xmin = list_chr[index][random - 1][0];
						ymin = list_chr[index][random - 1][1];
						zmin = list_chr[index][random - 1][2];
						double dist1 = sqrt((i-xmin)*(i-xmin)+(j-ymin)*(j-ymin)+(k-zmin)*(k-zmin));
						if(dist1 > d1){
							continue;
						}
					}else{
						xmax = list_chr[index][random + 1][0];
						ymax = list_chr[index][random + 1][1];
						zmax = list_chr[index][random + 1][2];
                                                xmin = list_chr[index][random - 1][0];
                                                ymin = list_chr[index][random - 1][1];
                                                zmin = list_chr[index][random - 1][2];
						double dist1 = sqrt((i-xmin)*(i-xmin)+(j-ymin)*(j-ymin)+(k-zmin)*(k-zmin));
                                                double dist2 = sqrt((i-xmax)*(i-xmax)+(j-ymax)*(j-ymax)+(k-zmax)*(k-zmax));
						if(dist1 > d1 || dist2 > d1){
							continue;
						}
					}
					struct node *temp;
					temp = new(struct node);
					temp -> x = i;
					temp -> y = j;
					temp -> z = k;
					candidates.push_back(*temp);
					//cout<<"candidate "<<i<<" "<<j<<" "<<k<<" added"<<endl;
					delete temp;		
				}
			}
		}
	}
	return candidates;

}
vector <node> Simulation::around_pos(int zero,int rand1,int length_seq,double d_1){
	vector <node> candidates_pos_1;
	//cout<<d1<<endl;
	//exit(1);
	double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;
	/*if(zero==1){
		x = list.get_x_value(rand1+1);
		y = list.get_y_value(rand1+1);
		z = list.get_z_value(rand1+1);
	}*/
	//exit(1);
	if(zero==2){
		x = list2.get_x_value(rand1 + 1);
                y = list2.get_y_value(rand1 + 1);
                z = list2.get_z_value(rand1 + 1);
	}
	//cout<<rand1<<endl;
	//cout<<x<<": "<<y<<": "<<z<<endl;
	//exit(1);
        for(double i = x - 1; i <= x + 1; i = i + 1){
		if(i < 0 || i >= length_seq){continue;}//change length;
                for(double j = y - 1; j <= y + 1;j = j + 1){
                        if(j < 0 || j >= length_seq){continue;}
			for(double k = z - 1; k <= z + 1;k = k + 1){
				if(k < 0 || k >= length_seq){ continue;}
                                if(x == i && y == j && z == k && list2.get_x2(rand1 + 1) != 1){
                                        continue;
                                }
                                else{
                                        /*struct node *temp;
                                        temp = new(struct node);
                                        temp -> x = i;
                                        temp -> y = j;
                                        temp -> z = k;*/
                                        if(rand1 + 1 > 1 && rand1 + 1 < list2.num_nodes){
						/*if(zero==1){
					                xmin = this->list.get_x_value(rand1);
					                xmax = this->list.get_x_value(rand1);
					                ymin = this->list.get_y_value(rand1);
					                ymax = this->list.get_y_value(rand1);
					                zmin = this->list.get_z_value(rand1);
					                zmax = this->list.get_z_value(rand1);
						}*/
						if(zero==2){
							xmin = this -> list2.get_x_value(rand1 + 1 - 1);
		                			xmax = this -> list2.get_x_value(rand1 + 1 + 1);
                					ymin = this -> list2.get_y_value(rand1 + 1 - 1);
        					        ymax = this -> list2.get_y_value(rand1 + 1 + 1);
              						zmin = this -> list2.get_z_value(rand1 + 1 - 1);
                					zmax = this -> list2.get_z_value(rand1 + 1 + 1);
						}                                                  
                                                double dist1 = sqrt((i - xmin) * (i - xmin) + (j - ymin) * (j - ymin) + (k - zmin) * (k - zmin));
                                                double dist2 = sqrt((i - xmax) * (i - xmax) + (j - ymax) * (j - ymax) + (k - zmax) * (k - zmax));                                                
                                               
						if(dist1 > (d_1 / 1) || dist2 > (d_1 / 1)){
                                                        continue;
                                                }
                                        }
                                        if(rand1  + 1 == 1){
						/*if(zero==1){
                                                        xmax = this->list.get_x_value(rand1+1);
                                                        ymax = this->list.get_y_value(rand1+1);
                                                        zmax = this->list.get_z_value(rand1+1);
						}*/
						if(zero == 2){
                                                        xmax = this -> list2.get_x_value(rand1 + 1 + 1);
                                                        ymax = this -> list2.get_y_value(rand1 + 1 + 1);
                                                        zmax = this -> list2.get_z_value(rand1 + 1 + 1);
						}
						double dist2 = sqrt((i - xmax)*(i - xmax)+(j - ymax)*(j - ymax)+(k - zmax)*(k - zmax));
                                                if(dist2 > (d_1 / 1)){
                                                        continue;
                                                }
                                        }
                                        if(rand1 + 1 == list2.num_nodes){
						if(zero == 2){
                                                        xmin = this -> list2.get_x_value(rand1 + 1 - 1);
                                                        ymin = this -> list2.get_y_value(rand1 + 1 - 1);
                                                        zmin = this -> list2.get_z_value(rand1 + 1 - 1);
						}
						/*if(zero==1){
                                                        xmin = this->list.get_x_value(rand1-1);
                                                        ymin = this->list.get_y_value(rand1-1);
                                                        zmin = this->list.get_z_value(rand1-1);
						}*/
                                                double dist1 = sqrt((i - xmin)*(i - xmin)+(j - ymin)*(j - ymin)+(k - zmin)*(k - zmin));
						if(dist1 > (d_1 / 1)){
                                                        continue;
                                                }
                                        }
					struct node *temp;
					temp = new(struct node);
					temp -> x = i;
					temp -> y = j;
					temp -> z = k;
                                        candidates_pos_1.push_back(*temp);
                                        delete temp;
                                }
                        }
                }
        }
	return candidates_pos_1;
}
bool Simulation::simulate_n(int rand0, int cube_unit){
        vector < vector <node> > coord_o = online_list;
        vector < vector <node> > coord1 = online_list;
        //vector < vector <double> > olist_o = olist;
	//cout<<"m is "<<m<<endl;
	//cout<<"rand num is "<<rand0<<endl;
        double loss1 = calculate_loss1(rand0);
	//cout<<"cost 1 is "<<loss1<<endl;
        vector <node> chr;
	/*for(int i = 0; i < coord1[rand0].size(); i++){
		cout<<coord1[rand0][i].x<<" "<<coord1[rand0][i].y<<" "<<coord1[rand0][i].z<<endl;
	}
	cout<<"!!!!!!!!!!!!!!"<<endl;*/
        chr = rotate_coord(coord1[rand0], cube_unit);
	coord1[rand0] = chr;
        vector <node> chr112 = change_position1(coord1, coord1[rand0], 4, biggest_length);//4 is cube numbers as a side length;
	coord1[rand0] = chr112;
	//exit(1);
        online_list = coord1;
	//this -> online_list = online_list;
        /*for(int i = 0; i < online_list[rand0].size(); i++){
                cout<<online_list[rand0][i].x<<" "<<online_list[rand0][i].y<<" "<<online_list[rand0][i].z<<" "<<online_list[rand0][i].x2<<endl;
        }
	exit(1);*/
        double loss2 = calculate_loss1(rand0);
        //cout<<"cost 2 is "<<loss2<<endl;
	//exit(1);
	double Z = loss2 - loss1;
	//cout<<"cost is "<<Z<<endl;
        if(Z<=0){
                En = loss2;
                return true;
        }
        else{
                double rand3 = 0;
                timeval t1;
                gettimeofday(&t1, NULL);
                srand(t1.tv_usec * t1.tv_sec);
                rand3 = rand()%(999999+1)/(float)(999999+1);
                if(rand3 < exp(-Z / temperature)){
                        En = loss2;
                        return true;
                }else{
                        En = loss1;
                        online_list = coord_o;
                        return false;
                }
        }
}
bool Simulation::simulate_lite(int index, double temperature1){
	timeval t1;
	gettimeofday(&t1, NULL);
	srand(t1.tv_usec * t1.tv_sec);
	int rand1 = 0;
	//cout<<"this chr has "<<list_chr[index].size()<<endl;
	rand1 = rand() % list_chr[index].size();
	//cout<<"pick nodes: "<<rand1<<endl;
	vector <node> candidates_pos;
	candidates_pos = arounding(index, rand1);
	//cout<<"cand num is "<<candidates_pos.size();
	//int rand2 = rand() % candidates_pos.size();
	if(candidates_pos.size() == 0){
		return false;
	}
	//cout<<"created candidates_pos"<<endl;
	int rand2 = rand() % candidates_pos.size();
	double cost1 = calculate_cost2(index, rand1);
	double oldx = list_chr[index][rand1][0];
	double oldy = list_chr[index][rand1][1];
	double oldz = list_chr[index][rand1][2];
	list_chr[index][rand1][0] = candidates_pos[rand2].x;
	list_chr[index][rand1][1] = candidates_pos[rand2].y;
	list_chr[index][rand1][2] = candidates_pos[rand2].z;
	double cost2 = calculate_cost2(index, rand1);
	double Z = cost2 - cost1;
	if(Z <= 0){
		En = cost2;
		return true;
	}else{
		double rand3 = rand() % (999999 + 1) / (float)(999999 + 1);
		if(rand3 < exp(-Z / temperature1)){
			En = cost2;
			return true;
		}else{
			En = cost1;
			list_chr[index][rand1][0] = oldx;
			list_chr[index][rand1][1] = oldy;
			list_chr[index][rand1][2] = oldz;
			return false; 
		}
	}	
		
}
void Simulation::decrease(int index){
	int failed_times = 0;
	double temperature1 = 10;
	//cout<<index<<endl;
	int len = list_chr[index].size();
	//cout<<index<<" "<<len<<endl;
	int a_time = len * 10;
	int f_time = len * 100;
	for(int i = 1; i <= 1000000; i++){
		//cout<<index<<": tem: "<<temperature1<<endl;
		int counter = 0;
		int accept_counter = 0;
		while(counter < f_time){
			counter++;
			if(simulate_lite(index, temperature1) == true){
				accept_counter ++;
			}
			if(accept_counter == a_time){
				failed_times = 0;
				break;
			}
			if(counter == f_time - 1){
				failed_times++;
			}
		}
		temperature1 = temperature1 * 0.9;
		if(failed_times == 3){
			break;
		} 
	}
}

             
///////////////////////////////////////////////////////////////////////////////
vector <int> get_id(int id, int length){
        vector <int> v;
	int a = 0;
        int b = 0;
        int c = 0;
        for(int i = 0; i < length; i++){
                for(int j =0; j < length; j++){
			for(int k = 0; k < length; k++){
				int m = length*length*i + length*j + k;
				if(m == id){
					a = i;
					b = j;
					c = k;
					v.push_back(a);
					v.push_back(b);
					v.push_back(c);
				} 
			}
                }
        }
	return v;
}

vector < vector <double> > Simulation::filter(int r_length){
	int e = length_seq/r_length + 1;
	vector <int> x;
        vector <vector <double> > v3;
	for(int i = 0; i < r_length*r_length*r_length; i++){
		int a = 0;
                int b = 0;
                int c = 0;
                for(int j = 0; j < r_length; j++){
                        for(int k = 0; k < r_length; k++){
                                for(int l = 0; l < r_length; l++){
                                        int m =r_length*r_length*j + r_length*k + l;
                                        if(m == i){
                                                a = j;
                                                b = k;
                                                c = l;
                                        }
                                }
                        }
                }
		vector <vector <double> > v2;
        	for(int j = 0; j < list.num_nodes; j++){
			vector <double> v1;
        		if(list.get_x_value(j+1) >= e*a && list.get_x_value(j+1) < e*a + e && list.get_y_value(j+1) >= e*b && list.get_y_value(j+1) < e*b + e && list.get_z_value(j+1) >= e*c && list.get_z_value(j+1) < e*c + e ){
				//cout<<i<<" "<<j<<" "<<list.get_x_value(j+1)<<" "<<list.get_y_value(j+1)<<" "<<list.get_z_value(j+1)<<endl;
				v1.push_back(i);
				v1.push_back(j);
				v1.push_back(list.get_x_value(j+1));
				v1.push_back(list.get_y_value(j+1));
				v1.push_back(list.get_z_value(j+1));
				v2.push_back(v1);
				v3.push_back(v1);
                	}
		}
		if(v2.size() != 0){
			x.push_back(i);
			cout<<i<<" cube: x range is from "<<a*e<<" to "<<a*e+(e-1)<<endl;
                	cout<<"        y range is from "<<b*e<<" to "<<b*e+(e-1)<<endl;
                	cout<<"        z range is from "<<c*e<<" to "<<c*e+(e-1)<<endl;
		}
        	
        }
        return v3;
	
}

bool Simulation::simulate_cube(double temperature1){
        timeval t1;
        gettimeofday(&t1, NULL);
        srand(t1.tv_usec * t1.tv_sec);
	//cout<<cubes.size()<<" cubes "<<endl;
	int rand1 = rand() % cubes.size();//get_x_value should be i+1
	int rand2 = rand() % cubes[rand1].size();
	//cout<<rand2<<" node piced "<<endl;
	//cout<<"picked node is "<<cubes[rand1][rand2][1]<<endl;
	double cost1 = calculate_loss(rand1, rand2);
	//cout<<"cost1 : "<<cost1<<endl;
	int e2 = (biggest_length * 4);
	//cout<<"start"<<endl;
	vector <node> cand1 = around_pos(2, cubes[rand1][rand2][1], e2, d1);
	//cand1 = candidates_pos;
        if(cand1.size() == 0){
                //cout<<"!!! no candidates, returen false"<<endl;
                return false;
        }
	int rand3 = rand() % cand1.size();
	double oldx = list2.get_x_value(cubes[rand1][rand2][1] + 1);
	double oldy = list2.get_y_value(cubes[rand1][rand2][1] + 1);
	double oldz = list2.get_z_value(cubes[rand1][rand2][1] + 1);
	list2.set_value(cubes[rand1][rand2][1] + 1, cand1[rand3].x, cand1[rand3].y, cand1[rand3].z);
	/*for(int i = 0; i < cand1.size(); i++){
		cout<<cand1[i].x<<" "<<cand1[i].y<<" "<<cand1[i].z<<endl;
	}*/
	//exit(1);
	//cout<<"ok"<<endl;
	double cost2 = calculate_loss(rand1, rand2);

	//cout<<"before "<<oldx<<" "<<oldy<<" "<<oldz<<"; new "<<cand1[rand3].x<<" "<<cand1[rand3].y<<" "<<cand1[rand3].z<<endl;
	//cout<<"cost2 : "<<cost2<<endl;
	//exit(1);
	double Z = cost2 - cost1;
	if(Z <= 0){
		En = cost2;
		//cout<<"cost is "<<En<<endl;
		return true;
	}else{
		double rand3 = rand() % (999999+1)/(float)(999999+1);
		if(rand3 < exp(-Z/temperature1)){
			En = cost2;
		//	cout<<"cost is "<<En<<endl;
			return true;
		}else{
			En = cost1;
		//	cout<<"cost is "<<En<<endl;
			list2.set_value(cubes[rand1][rand2][1] + 1, oldx, oldy, oldz);
			return false;
		}
	}
}

void Simulation::simulated_annealing(double temperature){
	//mu.lock();
	vector <double> ls;
	/*if(nodes.size() == 1){
		return;
	}*/
	//Ls.clear();
	//En.clear();
	//ls.push_back(index);
	int time = list2.num_nodes * 10;
	int time2 = time / 10;
        int failed_times = 0;
	cout<<time<<endl;
	cout<<time2<<endl;
	double temperature1 = temperature;
        for(int i =  1; i <= 1000000; i++){
		cout<<"temperature: "<<temperature1<<endl;
                /*if(temperature1 < 1){
                        break;
                }*/
		/*if(i % 100 == 0){
			ls.push_back(En);	
		}*/
                int counter = 0;
                int accept_counter = 0;
                while(counter < time){
                        counter++;
			//cout<<En<<endl;
			//cout<<" this is "<<counter<<" th time step"<<endl;
                        if(simulate_cube(temperature1) == true){
                                accept_counter++;
                        }
                        if(accept_counter == time2){
                               failed_times = 0;
                               break;
                        }
                        if(counter == time - 1){
                                failed_times++;
                        }
			ls.push_back(En);
			//cout<<"En is "<<En<<" "<<endl;
			//sleep(5);
                }
		//cout<<"!!! is "<<En<<endl;
		//ls.push_back(En);
		//time_t my_time9 = time(NULL);
                //cout<<"@@@@@@@@time@@@@@@@:"<<ctime(&my_time9)<<endl;
                temperature1 = temperature1 * 0.65;
		if(temperature1 < 0.1){
			break;
		}
                if(failed_times == 3){
                        break;
                }
        }
                        
	Ls.push_back(ls);
	//mu.unlock();

}
void someJob(int x){
	int z = x*x;
	cout<<z<<endl;
}
void Simulation::cal_cube(int cube_num, int rand1, int rand2, int zero){
	double cost1 = 0.0;
        double cost1_1 = 0.0;
        double cost1_2 = 0.0;
        double cost1_3 = 0.0;
    if(zero == 1){
	for(int i = 0; i < cubes[cube_num].size(); i++){
		for(int j = 0; j < list2.num_nodes; j++){
			double x_i = list2.get_x_value(cubes[cube_num][i][1] + 1);
                        double y_i = list2.get_y_value(cubes[cube_num][i][1] + 1);
                        double z_i = list2.get_z_value(cubes[cube_num][i][1] + 1);
                        double x_j = list2.get_x_value(j + 1);
                        double y_j = list2.get_y_value(j + 1);
                        double z_j = list2.get_z_value(j + 1);
			//cost_out<<"cube is "<<cube_num<<"; node i is "<<cubes1[cube_num][i][1]<<": x y z "<<x_i<<" "<<y_i<<" "<<z_i<<"; node j is "<<j<<": x1 y1 z1 "<<x_j<<" "<<y_j<<" "<<z_j;
			double dist = sqrt( (x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j) + (z_i - z_j) * (z_i - z_j) );
			//cost_out<<"; dist is "<<dist<<"; hic is "<<real_dist2[cubes1[cube_num][i][1]][cubes1[cube_num][j][1]]<<"; delta is "<<delta2[cubes1[cube_num][i][1]][cubes1[cube_num][j][1]]<<"; cost is "<<cost1_1<<" + "<<cost1_2<<" + "<<cost1_3<<endl;
			
	                /*if(real_dist[cubes[cube_num][i][1]][j] == 1 || delta2[cubes[cube_num][i][1]][j] == 1){
                                cost1_1 = cost1_1 + ((dist-delta0)*(dist-delta0))/(delta0*delta0);
                        }
                        if(real_dist[cubes[cube_num][i][1]][j] == 0 && delta2[cubes[cube_num][i][1]][j] < 1 && delta2[cubes[cube_num][i][1]][j] > theta1){
                                double t1 = delta2[cubes[cube_num][i][1]][j];
                                t1 = pow(t1,1.0/3);
                                t1 = delta0/t1;
                                cost1_2 = cost1_2 + beta*(1- exp(-(dist-t1)*(dist-t1)/mu1));
                        }
                        if(real_dist[cubes[cube_num][i][1]][j]==0 && delta2[cubes[cube_num][i][1]][j] <= theta1){
                                double t2 = 0;
                                t2 = pow(theta1,1.0/3);
                                t2 = delta0/t2;
                                cost1_3 = cost1_3 + tau*(1- 1/(1+exp(-(dist-(t2-rho)))/phi)) ;
                        }*/
			//cost_out<<"cube "<<cube_num<<" cost is "<<cost1_1<<" + "<<cost1_2<<" + "<<cost1_3<<"; cube is "<<cube_num<<"; node i is "<<cubes1[cube_num][i][1]<<": x y z "<<x_i<<" "<<y_i<<" "<<z_i<<"; node j is "<<cubes1[cube_num][j][1]<<": x1 y1 z1 "<<x_j<<" "<<y_j<<" "<<z_j<<endl;
			//cost_out.close();

		}
	}
    }else{
	double x_i = list2.get_x_value(cubes[rand1][rand2][1] + 1);
        double y_i = list2.get_y_value(cubes[rand1][rand2][1] + 1);
        double z_i = list2.get_z_value(cubes[rand1][rand2][1] + 1);
	//double x22 = list2.get_x2(cubes[rand1][rand2][1] + 1);
	int x1 = cubes[rand1][rand2][1];
	//cout<<"node is "<<cubes[rand1][rand2][1]<<endl;
	//cout<<rand1<<" "<<cube_num<<endl;
      
            for(int i = 0; i < cubes[cube_num].size(); i++){
		//cout<<"cmp node "<<cubes[cube_num][i][1]<<endl;
		  int y1 = cubes[cube_num][i][1];
                //for(int j = i + 1; j < cubes[cube_num].size(); j++){
                  if(x1 != y1){
                        double x_j = list2.get_x_value(cubes[cube_num][i][1] + 1);
                        double y_j = list2.get_y_value(cubes[cube_num][i][1] + 1);
                        double z_j = list2.get_z_value(cubes[cube_num][i][1] + 1);
			//double x_22 = list2.get_x2(cubes[cube_num][i][1] + 1);
		    //if(x_22 != 1){
                        double dist = sqrt((x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j) + (z_i - z_j) * (z_i - z_j) );
			//cout<<"dist is "
			//cost_out<<"cube is "<<cube_num<<"; node i is "<<cubes1[rand1][rand2][1]<<": x y z "<<x_i<<" "<<y_i<<" "<<z_i<<"; node j is "<<cubes1[cube_num][i][1]<<": x1 y1 z1 "<<x_j<<" "<<y_j<<" "<<z_j<<"; dist is "<<dist<<"; hic is "<<real_dist3[cubes1[rand1][rand2][1]][cubes1[cube_num][i][1]]<<"; delta is "<<delta3[cubes1[rand1][rand2][1]][cubes1[cube_num][i][1]]<<"; cost is "<<cost1_1<<" + "<<cost1_2<<" + "<<cost1_3<<endl;
                        if( delta2[x1][y1] < 1 && delta2[x1][y1] >= 0.95){
                                cost1_1 = cost1_1 + ((dist-delta0)*(dist-delta0))/(delta0*delta0);
                        }
                        if(delta2[x1][y1] < 0.95 && delta2[x1][y1] >= 0.7){
                                double t1 = delta2[x1][y1];
                                t1 = pow(t1,1.0/3);
                                t1 = delta0/t1;
                                cost1_2 = cost1_2 + beta*(1- exp(-(dist-t1)*(dist-t1)/mu1));
                        }
                        if(delta2[x1][y1] < 0.7){
                                double t2 = 0;
                                t2 = pow(0.7,1.0/3);
                                t2 = delta0/t2;
                                cost1_3 = cost1_3 + tau * (1 - (1 / (1 + exp( - (dist - (t2 - rho))) / phi)));
                        }
			//cost_out<<"cube "<<cube_num<<" cost is "<<cost1_1<<" + "<<cost1_2<<" + "<<cost1_3<<endl;
		    //}
                }
        }
      
    }
	cost1 = cost1_1 + cost1_2 + cost1_3;
	//cost_out<<"all cost is "<<cost1<<endl;
	//cost_out.close();
	//cout<<cube_num<<"!!!! cube is "<<cost1<<endl;
	cube_cost[cube_num] = cost1;
}
void Simulation::chr_sum(int rand1, int index, char *file_name){
	double cost = 0.0;
	double cost2 = 0.0;
	double cost3 = 0.0;
	double cost1 = 0.0;
	double d1sprt = d1 * d1;
        //ofstream cost_out;
        //cost_out.open(file_name, fstream::app);
        //cout<<"!~@#$"<<rand1<<endl;
    if(rand1 != index){
	for(int i = 0; i < online_list[rand1].size(); i++){
                double x1 = online_list[rand1][i].x;
                double y1 = online_list[rand1][i].y;
                double z1 = online_list[rand1][i].z;
	    	double x2 = online_list[rand1][i].x2;
		//cout<<"get in i is "<<i<<endl;
		//if(x2 != 1){
		    for(int j = 0; j < online_list[index].size(); j++){
			double x = online_list[index][j].x;
	                double y = online_list[index][j].y;
        	        double z = online_list[index][j].z;
			double x22 = online_list[index][j].x2;
		      //if(x22 != 1){
			double dist = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1) + (z - z1)*(z - z1));
			//cost_out<<"pairs "<<i<<" "<<j<<": "<<x1<<" "<<y1<<" "<<z1<<" "<<x<<" "<<y<<" "<<z<<endl;
			int i1 = list_start_end[rand1][0] + i;
			int j1 = list_start_end[index][0] + j;
			//cout<<i1<<" "<<j1<<" "<<dist<<endl;
			if(delta2[i1][j1] < 1 && delta2[i1][j1]>=0.95){
				//cost_out<<"pairs "<<i<<" "<<j<<" dist is "<<dist<<" delta0 is "<<delta0<<endl;
				cost1 = cost1 + ((dist-delta0)*(dist-delta0))/(delta0*delta0);
				//cost_out<<"cost is "<<cost1<<endl;
			}
			if(delta2[i1][j1] < 0.95 && delta2[i1][j1] >= 0.7){
				double t1 = delta2[i1][j1];
				t1 = pow(t1, 1.0/3);
				t1 = delta0/t1;
				cost2 = cost2 + beta*(1 - exp(- (dist - t1) * (dist - t1) / mu1));
			}
			if(delta2[i1][j1] < 0.7){
				double t2 = 0;
				t2 = pow(0.7, 1.0/3);
				t2 = delta0/t2;
				cost3 = cost3 + tau*(1 - 1/(1 + exp(-(dist - (t2 - rho)))/phi));
			}
			//cost_out<<"cost is "<<cost1<<" "<<cost2<<" "<<cost3<<endl;		
		      //}
			//}	
		    }
		//}
	    //}
	}
    }
	cost = cost1 + cost2 + cost3;
	rearrenge_cost[index] = cost;
	//cout<<"OK"<<endl;
	//return cost;
}
void Simulation::call_thread4(double rand1, double rand2){
        const size_t THREAD_COUNT = cubes.size();
	//cout<<cubes.size()<<" in call thread4 size is "<<endl;	
        std::vector<std::thread> threadPool;
        std::mutex mtx;
        for(int i = 0; i < cubes.size(); i++){
            //if(i != olist2[rand1][0] - 1){//from here;
                {
                 	//cout<<i<<" in call_thread4 "<<endl;
		        std::lock_guard<std::mutex> lock(mtx);
                        threadPool.emplace_back(std::thread([this, &mtx, &threadPool, i, rand1, rand2](){
                                //chr_sum(olist2, o_l, index_l, x_l, rand1, i);
                                cal_cube(i, rand1, rand2, 0);
                                std::lock_guard<std::mutex> lock(mtx);
                                threadPool.erase(
                                        std::find_if(threadPool.begin(), threadPool.end(),
                                        [](std::thread& x)
                                        {
                                                if(x.get_id() == std::this_thread::get_id()){
                                                        x.detach();
                                                        return true;
                                                }
                                                return false;
                                        })
                                );
                        }));
                }
                for(;;){
                        std::this_thread::yield();
                        std::lock_guard<std::mutex> lock(mtx);
                        if(threadPool.size() < THREAD_COUNT){
                                break;
                        }
                }
            //}
        }
        for(;;){
                std::this_thread::yield();
                std::lock_guard<std::mutex> lock(mtx);
                if(threadPool.empty()){
                        break;
                }
        }
}

void Simulation::call_thread3(int rand1){
        const size_t THREAD_COUNT = online_list.size();//need to change soon; 
	std::vector<std::thread> threadPool;
        std::mutex mtx;
        for(int i = 0; i < online_list.size(); i++){
	    //if(i !=  rand1){//from here;
                {
                        std::lock_guard<std::mutex> lock(mtx);
                        threadPool.emplace_back(std::thread([this, &mtx, &threadPool, rand1, i](){
				string cube_num = to_string(rand1);
				string cube_num2 = to_string(i);
				string file_o_name = "./cube_print/part" + cube_num + "_" + cube_num2;
				char *file_out = const_cast<char*>(file_o_name.c_str());
				chr_sum(rand1, i, file_out); 
                        	std::lock_guard<std::mutex> lock(mtx);
                                threadPool.erase(
                                        std::find_if(threadPool.begin(), threadPool.end(),
                                        [](std::thread& x)
                                        {
                                                if(x.get_id() == std::this_thread::get_id()){
                                                        x.detach();
                                                        return true;
                                                }
                                                return false;
                                        })
                                );
                        }));
                }
                for(;;){
                        std::this_thread::yield();
                        std::lock_guard<std::mutex> lock(mtx);
                        if(threadPool.size() < THREAD_COUNT){
                                break;
                        }
                }
            //}
        }
        for(;;){
                std::this_thread::yield();
                std::lock_guard<std::mutex> lock(mtx);
                if(threadPool.empty()){
                        break;
                }
        }
}

void Simulation::call_thread2(int chr_size){
        const size_t THREAD_COUNT = chr_size;
        std::vector<std::thread> threadPool;
        std::mutex mtx;
        for(int i = 0; i < chr_size; i++){
                {
                        std::lock_guard<std::mutex> lock(mtx);
                        threadPool.emplace_back(std::thread([this, &mtx, &threadPool, i](){
				decrease(i);               
				std::lock_guard<std::mutex> lock(mtx);
                                threadPool.erase(
                                        std::find_if(threadPool.begin(), threadPool.end(),
                                        [](std::thread& x)
                                        {
                                                if(x.get_id() == std::this_thread::get_id()){
                                                        x.detach();
                                                        return true;
                                                }
                                                return false;
                                        })
                                );
                        }));
                }
                for(;;){
                        std::this_thread::yield();
                        std::lock_guard<std::mutex> lock(mtx);
                        if(threadPool.size() < THREAD_COUNT){
                                break;
                        }
                }
        }
        for(;;){
                std::this_thread::yield();
                std::lock_guard<std::mutex> lock(mtx);
                if(threadPool.empty()){
                        break;
                }
        }
}

void Simulation::call_thread(vector <int> x1, vector <vector <int> > nodes_in1, vector <vector <double> > b1, int opz){
	Ls.clear();
    if(opz == 0){
	const size_t THREAD_COUNT = 10;
	std::vector<std::thread> threadPool;
	std::mutex mtx;
	for(int i = 0; i < x1.size(); i++){
		{
			int index1 = x1[i];
			vector <int> nodes_in2 = nodes_in1[i];
			cout<<nodes_in2.size()<<endl;
			vector <double> b2 = b1[i];
			std::lock_guard<std::mutex> lock(mtx);
			threadPool.emplace_back(std::thread([this, &mtx, &threadPool, index1, nodes_in2, b2](){
				//simulated_annealing();
				std::lock_guard<std::mutex> lock(mtx);
				threadPool.erase(
					std::find_if(threadPool.begin(), threadPool.end(),
					[](std::thread& x)
					{
						if(x.get_id() == std::this_thread::get_id()){
							x.detach();
							return true;
						}
						return false;
					})
				);
			}));
		}
		for(;;){
			std::this_thread::yield();
			std::lock_guard<std::mutex> lock(mtx);
			if(threadPool.size() < THREAD_COUNT){
				break;
			}
		}
	}
	for(;;){
		std::this_thread::yield();
		std::lock_guard<std::mutex> lock(mtx);
		if(threadPool.empty()){
			break;
		}
	}
   }else{
	//simulated_annealing(x1[0], nodes_in1[0], b1[0], 1);
   }
}

void Simulation::no_move(int r_length, double temperature){
	/*int e = ( biggest_length * 3)/r_length + 1;
	if(r_length == 1){
		e = biggest_length * 3;
	}*/
	//cout<<"e: length of the big cube is "<<e<<endl;
        vector <int> x;
	vector <int> x_b;
        vector <vector <double> > v3;
        vector <vector <double> > b1;
	vector <vector <double> > b1_b;
	vector <vector <int> > nodes_in1;
	vector <vector <int> > nodes_in1_b;
	vector < vector <vector <int> > > cube_list;
	//Double_list list2;
	//this -> list2 = list2;
        /*for(int i = 0; i < online_list.size(); i++){
                for(int k = 0; k < olist.size(); k++){
                        if(i == olist[k][0]){
                                for(int j = 0; j < online_list[k].size(); j++){
                                        list2.append_node(online_list[k][j].x, online_list[k][j].y, online_list[k][j].z);
                                }
                        }
                }
        }*/
	int len2 = (list2.num_nodes / np) + 1;
	int len3 = np;
	//cout<<len2<<" "<<len3<<"len2 len3"<<endl;
	for(int i = 0; i < np; i++){//generate a null list;
		vector <vector <int> > l1;
		cube_list.push_back(l1);
	}
	cout<<list2.num_nodes<<" nodes in list"<<endl;
	for(int i = 0; i < list2.num_nodes; i++){
		//int x1 = int( list2.get_x_value(i + 1) / e) + 1;
		//int y1 = int( list2.get_y_value(i + 1) / e) + 1;
		//int z1 = int( list2.get_z_value(i + 1) / e) + 1;
		//int in = (x1 - 1) * r_length * r_length + (y1 - 1) * r_length + (z1 - 1);
		int in = (i / len2); 
		//cout<<"i and in :"<<i<<" "<<in<<endl;
		vector <int> l1;
		l1.push_back(in);
		l1.push_back(i);
		l1.push_back(list2.get_x_value(i + 1));
		l1.push_back(list2.get_y_value(i + 1));
		l1.push_back(list2.get_z_value(i + 1));
		cout<<i<<endl;
		//l1.push_back(x1);
		//l1.push_back(y1);
		//l1.push_back(z1);
		cube_list[in].push_back(l1);
		//cout<< int( list2.get_x_value(i + 1) / e) + 1 <<" "<< int( list2.get_y_value(i + 1) / e) + 1<<" "<< int( list2.get_z_value(i + 1) / e) + 1<<endl;
	}
	//cout<<cube_list.size()<<endl;
	//exit(1);
	cout<<"original move"<<endl;
	int cube_count = 0;
	for(int i = 0; i < len3; i++){
		if(cube_list[i].size() != 0){
			vector <vector <int> > cube_lite;
			for(int j = 0; j < cube_list[i].size(); j++){
				vector <int> cube_lite_node;
				cube_lite_node.push_back(cube_count);//cube number
				cube_lite_node.push_back(cube_list[i][j][1]);//nodes number
				cube_lite.push_back(cube_lite_node);
			}
			cube_count++;
			cubes.push_back(cube_lite);	
		}
	}
	cube_cost.clear();
	for(int i = 0; i < cubes.size(); i++){
		cube_cost.push_back(0);
	}
	simulated_annealing(temperature);
}
void Simulation::move_boundary(int r_length){
	int e = (biggest_length * 3) / r_length + 1;
	int e2 = (e / 2) + 1;
	vector <int> x;
	vector <int> x_b;
	vector <int> x_c;
	vector <vector <double> > b1;
	vector <vector <double> > b1_b;
	vector <vector <double> > b1_c;
	vector <vector <vector <int> > > cube_list;
	vector <vector <vector <int> > > cube_list2;
	vector <vector <vector <int> > > cube_list3;
	vector <vector <int> > nodes_in1;
	vector <vector <int> > nodes_in1_b;
	vector <vector <int> > nodes_in1_c;
        for(int i = 0; i < r_length * r_length * (r_length + 1); i++){
                vector <vector <int> > l1;
                cube_list.push_back(l1);
	}
        for(int i = 0; i < list2.num_nodes; i++){
                int x1 = int( list2.get_x_value(i + 1) / e2) + 1;
		if(x1 != 1){
			if(x1 == 2 * r_length){
				x1 = r_length + 1;
			}else{
				x1 = int(x1 / 2) + 1;
			}
		}
                int y1 = int( list2.get_y_value(i + 1) / e) + 1;
                int z1 = int( list2.get_z_value(i + 1) / e) + 1;
                int in = (x1 - 1) * (r_length + 1) * r_length + (y1 - 1) * r_length + (z1 - 1);
                vector <int> l1;
                l1.push_back(in);
                l1.push_back(i);
                l1.push_back(list2.get_x_value(i + 1));
                l1.push_back(list2.get_y_value(i + 1));
                l1.push_back(list2.get_z_value(i + 1));
                l1.push_back(x1);
                l1.push_back(y1);
                l1.push_back(z1);
                cube_list[in].push_back(l1);
        }
	cout<<"x move "<<endl;
        for(int i = 0; i < r_length * r_length * (r_length + 1); i++){
                if(cube_list[i].size() != 0){
                        vector <int> nodes;
                        vector <double> boundary;
                        //cout<<"cube "<<cube_list[i][0][0]<<" has "<<cube_list[i].size()<<" nodes"<<endl;
                        x.push_back(cube_list[i][0][0]);
			if(cube_list[i][0][5] == 1){
				boundary.push_back(0);
				boundary.push_back(e2);
			}else if(cube_list[i][0][5] == r_length + 1){
				boundary.push_back(e2 + (cube_list[i][0][5] - 1) * e);
				boundary.push_back(2 * e2 + (cube_list[i][0][5] - 1) * e);
			}else{
				boundary.push_back(e2 + (cube_list[i][0][5] - 1) * e);
				boundary.push_back(e2 + (cube_list[i][0][5] ) * e);
			}
                        boundary.push_back((cube_list[i][0][6] - 1) * e);
                        boundary.push_back( cube_list[i][0][6] * e);
                        boundary.push_back((cube_list[i][0][7] - 1) * e);
                        boundary.push_back( cube_list[i][0][7] * e);
			b1.push_back(boundary);
                        for(int j = 0; j < cube_list[i].size(); j++){
                                nodes.push_back(cube_list[i][j][1]);
                        }
                        nodes_in1.push_back(nodes);
                }
        }
        /*for(int i = 0; i < nodes_in1.size(); i++){
                for(int j = 0; j < nodes_in1[i].size(); j++){
                        cout<<b1[i][0]<<" "<<b1[i][1]<<"> <"<<b1[i][2]<<" "<<b1[i][3]<<"> <"<<b1[i][4]<<" "<<b1[i][5]<<"> "<<list2.get_x_value(nodes_in1[i][j] + 1)<<" "<<list2.get_y_value(nodes_in1[i][j] + 1)<<" "<<list2.get_z_value(nodes_in1[i][j] + 1)<<endl;
                }

        }
	exit(1);*/
        cout<<"start x multithreading#####"<<endl;
        call_thread(x, nodes_in1, b1, 0);

        string s1212= "2_11_xstep_e5";//store the delata matrix
        char *p52= const_cast<char*>(s1212.c_str());
        ofstream file1572;
        file1572.open(p52);
        for(int i = 0; i < Ls.size(); i++){
                for(int j = 0; j < Ls[i].size();j++){
                        file1572<<Ls[i][j]<<" ";
                }
                file1572<<endl;
        }
        file1572.close();
        cout<<"finish x multithreading###"<<endl;
		
        for(int i = 0; i < r_length * r_length * (r_length + 1); i++){
                vector <vector <int> > l1;
                cube_list2.push_back(l1);
        }
	cout<<cube_list2.size();
        for(int i = 0; i < list2.num_nodes; i++){
                int x1 = int( list2.get_x_value(i + 1) / e) + 1;
                int y1 = int( list2.get_y_value(i + 1) / e2) + 1;
		if(y1 != 1){
			if(y1 == 2 * r_length){
				y1 = r_length + 1;
			}else{
                        	y1 = int(y1 / 2) + 1;
			}
                }
                int z1 = int( list2.get_z_value(i + 1) / e) + 1;
		int in = 0;
		if(y1 != r_length + 1){
                	in = (x1 - 1) * (r_length) * r_length + (y1 - 1) * r_length + (z1 - 1);
                }else{
			in = r_length * r_length * r_length + (x1 - 1) * r_length + (z1 - 1);
		}
		//cout<<"cube "<<in<<" is "<<x1<<" "<<y1<<" "<<z1<<endl;
		vector <int> l1;
                l1.push_back(in);
                l1.push_back(i);
                l1.push_back(list2.get_x_value(i + 1));
                l1.push_back(list2.get_y_value(i + 1));
                l1.push_back(list2.get_z_value(i + 1));
                l1.push_back(x1);
                l1.push_back(y1);
                l1.push_back(z1);
                cube_list2[in].push_back(l1);
        }
	//exit(1);
        cout<<"y move "<<endl;
        for(int i = 0; i < r_length * r_length * (r_length + 1); i++){
                if(cube_list2[i].size() != 0){
                        vector <int> nodes;
                        vector <double> boundary;
                        //cout<<"cube "<<cube_list2[i][0][0]<<" has "<<cube_list2[i].size()<<" nodes"<<endl;
                        x_b.push_back(cube_list2[i][0][0]);
                        boundary.push_back((cube_list2[i][0][5] - 1) * e);
                        boundary.push_back( cube_list2[i][0][5] * e);
			cout<<"y1 "<<cube_list2[i][0][6]<<endl;
                        if(cube_list2[i][0][6] == 1){
                                boundary.push_back(0);
                                boundary.push_back(e2);
                        }else if(cube_list2[i][0][6] == r_length + 1){
                                boundary.push_back(e2 + (cube_list2[i][0][6] - 1) * e);
                                boundary.push_back(2 * e2 + (cube_list2[i][0][6] - 1) * e);
                        }else{
                                boundary.push_back(e2 + (cube_list2[i][0][6] - 1) * e);
                                boundary.push_back(e2 + (cube_list2[i][0][6] ) * e);
                        }
                        boundary.push_back((cube_list2[i][0][7] - 1) * e);
                        boundary.push_back( cube_list2[i][0][7] * e);
                        b1_b.push_back(boundary);
                        for(int j = 0; j < cube_list2[i].size(); j++){
                                nodes.push_back(cube_list2[i][j][1]);
                        }
                        nodes_in1_b.push_back(nodes);
                }
        }
	/*for(int i = 0; i < nodes_in1_b.size(); i++){
		for(int j = 0; j < nodes_in1_b[i].size(); j++){
			cout<<b1_b[i][0]<<" "<<b1_b[i][1]<<"> <"<<b1_b[i][2]<<" "<<b1_b[i][3]<<"> <"<<b1_b[i][4]<<" "<<b1_b[i][5]<<"> "<<list2.get_x_value(nodes_in1_b[i][j] + 1)<<" "<<list2.get_y_value(nodes_in1_b[i][j] + 1)<<" "<<list2.get_z_value(nodes_in1_b[i][j] + 1)<<endl;
		}
		
	}*/
	//exit(1);
        cout<<"start y multithreading#####"<<endl;
        call_thread(x_b, nodes_in1_b, b1_b, 0);
        string s12121= "2_11_ystep_e5";//store the delata matrix
        char *p521= const_cast<char*>(s12121.c_str());
        ofstream file15721;
        file15721.open(p521);
        for(int i = 0; i < Ls.size(); i++){
                for(int j = 0; j < Ls[i].size();j++){
                        file15721<<Ls[i][j]<<" ";
                }
                file15721<<endl;
        }
        file15721.close();
        cout<<"finish y multithreading###"<<endl;
        for(int i = 0; i < r_length * r_length * (r_length + 1); i++){
                vector <vector <int> > l1;
                cube_list3.push_back(l1);
        }
        for(int i = 0; i < list2.num_nodes; i++){
                int x1 = int( list2.get_x_value(i + 1) / e) + 1;
                int y1 = int( list2.get_y_value(i + 1) / e) + 1;
                int z1 = int( list2.get_z_value(i + 1) / e2) + 1;
		if(z1 != 1){
			if(z1 == 2 * r_length){
				z1 = r_length + 1;
			}else{
				z1 = int(z1 / 2) + 1;
			}
		}
		int in = 0;
                if(z1 != r_length + 1){
                        in = (x1 - 1) * (r_length) * r_length + (y1 - 1) * r_length + (z1 - 1);
                }else{
                        in = r_length * r_length * r_length + (x1 - 1) * r_length + (y1 - 1);
                }
                vector <int> l1;
                l1.push_back(in);
                l1.push_back(i);
                l1.push_back(list2.get_x_value(i + 1));
                l1.push_back(list2.get_y_value(i + 1));
                l1.push_back(list2.get_z_value(i + 1));
                l1.push_back(x1);
                l1.push_back(y1);
                l1.push_back(z1);
                cube_list3[in].push_back(l1);
        }
        cout<<"z move "<<endl;
        for(int i = 0; i < r_length * r_length * (r_length + 1); i++){
                if(cube_list3[i].size() != 0){
                        vector <int> nodes;
                        vector <double> boundary;
			//cout<<"cube "<<cube_list3[i][0][0]<<" has "<<cube_list3[i].size()<<" nodes"<<endl;
                        x_c.push_back(cube_list3[i][0][0]);
                        boundary.push_back((cube_list3[i][0][5] - 1) * e);
                        boundary.push_back( cube_list3[i][0][5] * e);
                        boundary.push_back((cube_list3[i][0][6] - 1) * e);
                        boundary.push_back( cube_list3[i][0][6] * e);
                        if(cube_list3[i][0][7] == 1){
                                boundary.push_back(0);
                                boundary.push_back(e2);
                        }else if(cube_list3[i][0][7] == r_length + 1){
                                boundary.push_back(e2 + (cube_list3[i][0][7] - 1) * e);
                                boundary.push_back(2 * e2 + (cube_list3[i][0][7] - 1) * e);
                        }else{
                                boundary.push_back(e2 + (cube_list3[i][0][7] - 1) * e);
                                boundary.push_back(e2 + (cube_list3[i][0][7] ) * e);
                        }
                        b1_c.push_back(boundary);
                        for(int j = 0; j < cube_list3[i].size(); j++){
                                nodes.push_back(cube_list3[i][j][1]);
                        }
                        nodes_in1_c.push_back(nodes);
                }
        }
        cout<<"start z multithreading#####"<<endl;
        call_thread(x_c, nodes_in1_c, b1_c, 0);
        string s12123= "2_11_zstep_e5";//store the delata matrix
        char *p522= const_cast<char*>(s12123.c_str());
        ofstream file15722;
        file15722.open(p522);
        for(int i = 0; i < Ls.size(); i++){
                for(int j = 0; j < Ls[i].size();j++){
                        file15722<<Ls[i][j]<<" ";
                }
                file15722<<endl;
        }
        file15722.close();
        cout<<"finish z multithreading###"<<endl;
 //       cout<<"after simulation"<<endl;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
}
/////////////////////////////////////////////////////////////////////////////////////////int main //////////////////////
int main(int argc, char* argv[]){
	int seq_length = 0;
	double temperature = 10.0;
	string dist_fi = "none";
	string dist_fi2 = "none";
	string impute_fi = "nonei";
	string or_fi = "noner";
	string out_f = "none1";
	string out_f2 = "none2";
	string out_f3 = "none3";
	string list_name = "mouse_len";
	double rho = 1.0;
	double phi = 0.1;
	double theta1 = 0.7;
	double beta = 1.0;
	double delta0 = 8.0;
	double mu1 = 20.0;
	double mu2 = 2.0;
	double tau = 1.0;
	int d0 = 5;
	int np = 40;
	int open = 0;
	double d1 = 8.0;
	int seq_length2 = 0;
	double res = 0.5;
	double temperature2 = 0.1;
	double e1 = 0.24;
	int th = 7;	
	for(int i = 1; i < argc - 1; i++){
                if(strcmp(argv[i],"-i")==0){
                        dist_fi = argv[i+1];
                }else if(strcmp(argv[i],"-i2")==0){
			dist_fi2 = argv[i+1];
		}else if(strcmp(argv[i],"-o")==0){
                        out_f = argv[i+1];
                }else if(strcmp(argv[i],"-len")==0){
			list_name = argv[i+1];
		}else if(strcmp(argv[i],"-o2")==0){
			out_f2 = argv[i+1];
		}else if(strcmp(argv[i],"-o3")==0){
			out_f3 = argv[i+1];
		}else if(strcmp(argv[i],"-l")==0){
                        string length (argv[i+1]);
                        seq_length = atoi(length.c_str());
                }else if(strcmp(argv[i],"-t")==0){
                        string temper (argv[i+1]);
                        temperature = atof(temper.c_str());
                }else if(strcmp(argv[i],"-rho")==0){
			string r (argv[i+1]);
			rho = atof(r.c_str()); 
		}else if(strcmp(argv[i],"-theta1")==0){
                        string r (argv[i+1]);
                        theta1 = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-delta0")==0){
                        string r (argv[i+1]);
                        delta0 = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-mu1")==0){
                        string r (argv[i+1]);
                        mu1 = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-beta")==0){
                        string r (argv[i+1]);
                        beta = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-tau")==0){
                        string r (argv[i+1]);
                        tau = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-phi")==0){
                        string r (argv[i+1]);
                        phi = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-d0")==0){//d0 now is the argument to modified the simulation times
                        string r (argv[i+1]);
                        d0 = atoi(r.c_str()); 
                }else if(strcmp(argv[i],"-np")==0){
			string r (argv[i+1]);
			np = atoi(r.c_str());
		}else if(strcmp(argv[i],"-th")==0){
			string r (argv[i+1]);
			th = atoi(r.c_str());
		}else if(strcmp(argv[i],"-mu2")==0){
                        string r (argv[i+1]);
                        mu2 = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-cost")==0){
			string r (argv[i+1]);
			open = atoi(r.c_str());
		}else if(strcmp(argv[i],"-d1")==0){
			string r (argv[i+1]);
			d1 = atof(r.c_str());
		}else if(strcmp(argv[i],"-res")==0){
			string r (argv[i+1]);
			res = atof(r.c_str());
		}else if(strcmp(argv[i],"-e")==0){
			string r (argv[i+1]);
			e1 = atof(r.c_str());
		}else if(strcmp(argv[i],"-t2")==0){//for the expand list simulation start from
			string r (argv[i+1]);
			temperature2 = atof(r.c_str());
		}else if(strcmp(argv[i],"-or")==0){
                        or_fi = argv[i+1];
                }else if(strcmp(argv[i],"-impute")==0){
                        impute_fi = argv[i+1];
                }
        }
	char *dist_file=const_cast<char*>(dist_fi.c_str());
	if(dist_fi=="none"){
                cout<<"you should insert a hic contacts file"<<endl;
                exit(1);
        }
	vector < vector <double> > real_dist;//to save the hic matrix
	vector < vector <double> > real_dist2;//to save the high resolution hic matrix
	vector < vector <double> > delta;//to save the delta matrix
	vector < vector <double> > delta2;//to save the high resolution delta matrix
	vector < vector < vector <double> > > hic_list;
	/////////////////////////////////////////////////////scl///////////////
        //string s2 = "t2";//t2 is imputed matrix, used for insteading delta
        char *fls1 = const_cast<char*>(impute_fi.c_str());
	if(impute_fi=="nonei"){
		cout<<"you should insert an imputed hic matrix"<<endl;
		exit(1);
	}
        ifstream flsr1;
        flsr1.open(fls1);
        string line3;
        vector <vector<double> > tmat3;
        double temp3;
        while(getline(flsr1, line3)){
                vector <double> row;
                istringstream istr1(line3);
                while(istr1 >> temp3){
                        row.push_back(temp3);
                }

                tmat3.push_back(row);
        }
        delta2 = tmat3;
	vector <string> chromosome_num;
	chromosome_num = read_chr_nums(dist_file);
        char *file_ln = const_cast<char*>(list_name.c_str());
        //char *fls1 = const_cast<char*>(impute_fi.c_str());
        if(list_name=="mouse_len"){
                cout<<"you should insert an chromosome length file"<<endl;
                exit(1);
        }
        ifstream fileln;
        fileln.open(file_ln);
        string line2_1;
        vector <vector <string> > name_list;
        string temp2;
        while(getline(fileln,line2_1)){
                vector <string> row;
                istringstream istr1(line2_1);
                while(istr1 >> temp2){
                        row.push_back(temp2);
                }
                name_list.push_back(row);
        }

	hic_list = read_w_hic(res, dist_file, chromosome_num, name_list);//test!!!!!!!!!!!
	cout<<"finish generating hic matrix list; list is "<<hic_list.size()<<endl;
        vector< vector <int> > null_line;
        for(int i = 0; i < hic_list.size(); i++){
                //cout<<"test "<<hic_list[i].size()<<endl;
                vector <int> null;
                for(int j = 0; j < hic_list[i].size(); j++){
                        int countj = 0;
                        for(int z = 0; z < hic_list[i][j].size(); z++){
                                if(hic_list[i][j][z] == 0){
                                        countj++;
                                }
                        }
                        if(countj == hic_list[i].size()){
                                null.push_back(j);
                                //cout<<j<<"null line"<<endl;
                        }
                }
                null_line.push_back(null);
        }
	vector < vector < vector <double> > > delta_list;//test!!!!!!!!!!!!!
	/*for(int i = 0; i < hic_list.size(); i++){
		vector < vector <double> > temp_delta = calculate_t2(hic_list[i].size(), hic_list[i], mu2, d0);
		cout<<temp_delta.size()<<endl;
		delta_list.push_back(temp_delta);
	}*/
	vector <vector <int> > sam_list;
	int ccu = 0;
	for(int i = 0; i < hic_list.size(); i++){
		string str = to_string(i+1);
		if(i == hic_list.size() - 1){
			str = "X";
		}
		vector <int> etp;
		etp.push_back(ccu);
		cout<<ccu;
		for(int j = 0; j < hic_list.size(); j++){
			if(str == chromosome_num[j]){
				//cout<<chromosome_num[j]<<" "<<hic_list[j].size()<<endl;
				ccu+=hic_list[j].size();
			}
		}
		etp.push_back(ccu);
		//cout<<" "<<ccu<<endl;
		sam_list.push_back(etp);
		//cout<<str<<endl;
	}
	for(int i = 0; i < hic_list.size(); i++){
		vector < vector <double> > del;
		int m = atoi(chromosome_num[i].c_str());
		m = m -1;
		if(chromosome_num[i] == "X"){
			m = hic_list.size() - 1;
		}
		//cout<<m<<" "<<sam_list[m][0]<<" "<<sam_list[m][1]<<endl;
		for(int j = sam_list[m][0]; j < sam_list[m][1]; j++){
			vector <double> de2;
			for(int k = sam_list[m][0]; k < sam_list[m][1]; k++){
				de2.push_back(delta2[j][k]);
			}
			del.push_back(de2);
		}
		delta_list.push_back(del);
	}
	delta2.clear();
	cout<<"finish creating delta matrix; size is "<<delta_list.size()<<endl;
	for(int i = 0; i < hic_list.size(); i++){
		cout<<chromosome_num[i]<<" "<<hic_list[i].size()<<" "<<delta_list[i].size()<<endl;
		//for(int j = 0; j < hic_list[i].size(); j++){
		//	cout<<hic_list[i][j].size()<<" "<<delta_list[i][j].size()<<endl;
		//}
	}
	//exit(1);
	vector < vector < vector <double> > > chr_list;
	int size_cube = 0;//set the size cube as number of beads times 5 to let the structure moves sufficiently
        int sim_times = 0;///testinggggggggggg 
	char *out_file=const_cast<char*>(out_f.c_str());
	char *out_file2=const_cast<char*>(out_f2.c_str());
	char *out_file3=const_cast<char*>(out_f3.c_str());
	vector < vector <double> > nchr_p;
	Simulation simulation(nchr_p,size_cube, seq_length, 2, sqrt(10), hic_list, temperature, delta_list,rho,theta1,mu1,mu2,tau,beta,delta0,d0,phi,d1,seq_length2,real_dist2,delta2,e1,np);
	//cout<<"Initialization of the 3D structure started."<<endl;
	//read_chr_list_length;
	/*for(int i = 0; i < name_list.size(); i++){
		cout<<name_list[i][0]<<" "<<name_list[i][1]<<endl;
	}*/
	for(int i = 0; i < chromosome_num.size(); i++){
		int length_c = hic_list[i].size();
		
		//cout<<chromosome_num[i]<<" "<<length_c<<endl;
		simulation.length_seq = length_c;
		simulation.size_cube = length_c * 5;
		bool initial_status = false;
		while(initial_status == false){
			initial_status = simulation.rand_initialize();
		}
		vector < vector <double> > chr_s_list;
		for(int j = 1; j <= simulation.list.num_nodes; j++){
			vector <double> chr_xyz;
			chr_xyz.push_back(simulation.list.get_x_value(j));
			chr_xyz.push_back(simulation.list.get_y_value(j));
			chr_xyz.push_back(simulation.list.get_z_value(j));
			chr_s_list.push_back(chr_xyz);
			//cout<<i+1<<" "<<chr_xyz[0]<<" "<<chr_xyz[1]<<" "<<chr_xyz[2]<<endl;
		}
		simulation.list_chr.push_back(chr_s_list);
		cout<<simulation.list.num_nodes<<endl;
		simulation.list.clear();
		cout<<simulation.list.num_nodes<<endl;
	}
        ofstream file13;
        file13.open(out_file3);
        /*for(int i = 0; i < simulation.list_chr.size(); i++){
                //file13<<"chr"<<chromosome_num[i]<<"##########"<<endl;
                for(int j = 0; j < simulation.list_chr[i].size(); j++){
                        for(int k = 0; k < simulation.list_chr[i][j].size(); k++){
                                file13<<simulation.list_chr[i][j][k]<<" ";
                        }
                        file13<<endl;
                }
        }*/
	simulation.call_thread2(chromosome_num.size());//running multiple threads scl;
	delta_list.clear();
	simulation.delta_list.clear();	
	hic_list.clear();
	simulation.hic_list.clear();
	for(int i = 0; i < simulation.list_chr.size(); i++){
		cout<<"chr"<<chromosome_num[i]<<"##########"<<endl;
		for(int j = 0; j < simulation.list_chr[i].size(); j++){
			for(int k = 0; k < simulation.list_chr[i][j].size(); k++){
				file13<<simulation.list_chr[i][j][k]<<" ";
			}
			file13<<endl;
		}
	}
	file13.close();
	time_t my_time0 = time(NULL);
        cout<<"@@@@@@@@time@@@@@@@:"<<ctime(&my_time0)<<endl;
	//exit(1);
	///////////////////////////////////////rearrange///////////////////////////////
	vector < vector < vector <double> > > coordinates;
	int file_size = simulation.list_chr.size();
	vector <int> files;
	vector <string> files2;
	for(int i = 0; i < file_size; i++){
		//cout<<i<<" "<<chromosome_num[i]<<endl;
		files2.push_back(chromosome_num[i]);
		//cout<<files2[i]<<endl;
		/*if(chromosome_num[i] == "1"){
			vector <string> pp1;
	                vector <int> null;

		}*/
	}
	//exit(1);
	vector < vector <double> > chr_p;
	vector < vector <string> > chr_pp;
	vector < vector <int> > cnull_line;
        //string s = "ranlist";
        char *fls = const_cast<char*>(or_fi.c_str());
        if(or_fi=="noner"){
                cout<<"you should insert a rank list of chromosome"<<endl;
                exit(1);
        }

        //char *fls = const_cast<char*>(s.c_str());
        ifstream flsr;
        flsr.open(fls);
        string line1;
       vector <  vector<string> > tmat;
        string temp1;
        while(getline(flsr, line1)){
                vector <string> row;
                istringstream istr1(line1);
                while(istr1 >> temp1){
                        row.push_back(temp1);
                }

                tmat.push_back(row);
        }

	for(int i = 0; i < hic_list.size(); i++){
		vector <double> p1;
		vector <string> pp1;
		vector <int> null;
		timeval t1;
		gettimeofday(&t1, NULL);
		srand(t1.tv_usec * t1.tv_sec);
		//int rand1 = rand() % 64;
		//cout<<rand1<<" th "<<endl;
		//double value = files[rand1];
		string value2;
		/*int pos = 0;
		if(i==0){
			value2 = "1";
			pos = 32;
		}else{
			value2 = "10";
			pos = rand() % 64;
			while(pos == 32){
				pos = rand() %64;
			}
		}*/
		//cout<<value2<<endl;
		///files.erase(files.begin() + rand1);
		//files2.erase(files2.begin() + rand1);
		//cout<<value<<endl;
		//cout<<value2<<endl;
		//p1.push_back(value);
		//p1.push_back(simulation.list_chr[value].size());
		//p1.push_back(i);
		value2 = tmat[i][0];
	        int rand2 = 0;
                for(int j = 0; j < chromosome_num.size(); j++){
                        if(value2 == chromosome_num[j]){
                                rand2 = j;
                        }
                }
		//value2 = tmat[i][0];
		pp1.push_back(value2);
		pp1.push_back(to_string(simulation.list_chr[rand2].size()));
		pp1.push_back(to_string(i));
		//null = null_line[rand2];
		for(int i = 0; i < null_line[rand2].size(); i++){
			//cout<<null_line[rand2][i]<<endl;
			null.push_back(null_line[rand2][i]);
		}
		cnull_line.push_back(null);
		//chr_p.push_back(p1);
		chr_pp.push_back(pp1);
		cout<<pp1[0]<<" "<<pp1[1]<<" "<<pp1[2]<<endl;
		//exit(1);
		coordinates.push_back(simulation.list_chr[rand2]);
	}
	//exit(1);//random move the
	int biggest_length = 0;
	for(int i = 0; i < hic_list.size(); i++){
                vector <double> x_a;
                vector <double> y_a;
                vector <double> z_a;
		cout<<coordinates[i].size()<<endl;
                for(int i1 = 0; i1 < coordinates[i].size(); i1++){
                        x_a.push_back(coordinates[i][i1][0]);
                        y_a.push_back(coordinates[i][i1][1]);
                        z_a.push_back(coordinates[i][i1][2]);
                }
                double blx = cal_big(x_a);
                double bly = cal_big(y_a);
                double blz = cal_big(z_a);
                if(biggest_length < blx){
                        biggest_length = blx;
                }
                if(biggest_length < bly){
                        biggest_length = bly;
                }
                if(biggest_length < blz){
                        biggest_length = blz;
                }

	}
	cout<<"biggest length is "<<biggest_length<<endl;
	//exit(1);  
	for(int i = 0; i < file_size; i++){
		vector <double> x_a;
		vector <double> y_a;
		vector <double> z_a;
		for(int i1 = 0; i1 < coordinates[i].size(); i1++){
			x_a.push_back(coordinates[i][i1][0]);
                        y_a.push_back(coordinates[i][i1][1]);
                        z_a.push_back(coordinates[i][i1][2]);
		}
		double mid_x = calculate_middle_point1(x_a);
		double mid_y = calculate_middle_point1(y_a);
		double mid_z = calculate_middle_point1(z_a);
                double move_x = int(biggest_length / 2) - mid_x;//10, 10, 10 is the middle point of a small cube
                double move_y = int(biggest_length / 2) - mid_y;
                double move_z = int(biggest_length / 2) - mid_z;
                for(int i1 = 0; i1 < coordinates[i].size(); i1++){
                        coordinates[i][i1][0] = coordinates[i][i1][0] + move_x;
                        coordinates[i][i1][1] = coordinates[i][i1][1] + move_y;
                        coordinates[i][i1][2] = coordinates[i][i1][2] + move_z;
			//cout<<coordinates[i][i1][0]<<" "<<coordinates[i][i1][1]<<" "<<coordinates[i][i1][2]<<endl;
                }
		
	}
        vector < vector <node> > all_list;
	int all_size = 0;
	int cube_length_num = 4;
	int cube_cln = cube_length_num * cube_length_num * cube_length_num;
	vector <int> list_index;
	for(int i = 0; i < cube_cln; i++){
		list_index.push_back(i);
	}
	vector <int> position_in_cube;
	for(int i = 0; i < file_size; i++){
		timeval t1;
		gettimeofday(&t1, NULL);
		srand(t1.tv_usec * t1.tv_sec);
		int rand1 = rand() % list_index.size();
		int value = list_index[rand1];
		list_index.erase(list_index.begin() + rand1);
		//position_in_cube.push_back(value);
		if(i==0){
                        position_in_cube.push_back(32);
                }else{
                        
                        int pos = rand() % 64;
                        while(pos == 32){
                                pos = rand() %64;
                        }
			position_in_cube.push_back(pos);
                }

	}
	for(int i = 0; i < file_size; i++){
                cout<<"get position "<<position_in_cube[i]<<endl;
        }
	cout<<file_size<<endl;
	//exit(1);
	simulation.biggest_length = biggest_length;
	for(int i = 0; i < file_size; i++){
                vector <node> one_list;
                all_size = all_size + coordinates[i].size();
                for(int i1 = 0; i1< coordinates[i].size(); i1++){

                        for(int j = 0; j < cube_length_num; j++){
                                for(int k = 0; k < cube_length_num; k++){
                                        for(int l = 0; l < cube_length_num; l++){
                                                if(position_in_cube[i] == cube_length_num * cube_length_num * j + cube_length_num * k + l){
                                                        
							coordinates[i][i1][0] = coordinates[i][i1][0] + j * biggest_length;
                                                        coordinates[i][i1][1] = coordinates[i][i1][1] + k * biggest_length;
                                                        coordinates[i][i1][2] = coordinates[i][i1][2] + l * biggest_length;

                                                        struct node*temp;
                                                        temp = new(struct node);
                                                        temp -> x = coordinates[i][i1][0];
                                                        temp -> y = coordinates[i][i1][1];
                                                        temp -> z = coordinates[i][i1][2];
                                                        one_list.push_back(*temp);
                                                        delete(temp);
                                                }
                                        }
                                }
                        }
                }
	//	cout<<one_list.size()<<endl;
		all_list.push_back(one_list);//testinggggggggg

	}
	for(int i = 0; i < cnull_line.size(); i++){
		cout<<"chr size "<<all_list[i].size()<<endl;
		for(int j = 0; j < cnull_line[i].size(); j++){
			all_list[i][cnull_line[i][j]].x2 = 1;
			//cout<<" null line "<<cnull_line[i][j]<<endl;
		}
	}
	/*for(int i = 0; i < all_list.size(); i++){
		for(int j = 0; j < all_list[i].size(); j++){
			cout<<all_list[i][j].x2<<endl;
		}
	}*/
	//simulation.cnull_line = cnull_line;
	/*cout<<all_list.size()<<endl;
	cout<<"################"<<endl;
	for(int i = 0; i < all_list.size(); i++){
		cout<<all_list[i].size()<<endl;
		for(int j = 0; j < all_list[i].size(); j++){
			cout<<all_list[i][j].x<<" "<<all_list[i][j].y<<" "<<all_list[i][j].z<<endl;
		}
	}*/
	//exit(1);
	cout<<"############"<<endl;
	//exit(1);*/
	for(int i = 0; i < chr_pp.size(); i++){
		for(int j = 0; j < chr_pp[i].size(); j++){
			cout<<chr_pp[i][j]<<" ";
		}
		cout<<endl;
	}
        /*string s2 = "t2";
        char *fls1 = const_cast<char*>(s2.c_str());
        ifstream flsr1;
        flsr1.open(fls1);
        string line3;
        vector <vector<double> > tmat3;
        double temp3;
        while(getline(flsr1, line3)){
                vector <double> row;
                istringstream istr1(line3);
                while(istr1 >> temp3){
                        row.push_back(temp3);
                }

                tmat3.push_back(row);
        }
        delta2 = tmat3;*/
        /*string s3 = "hic_896";
        char *fls2 = const_cast<char*>(s3.c_str());
        ifstream flsr2;
        flsr2.open(fls2);
        string line4;
        vector <vector<double> > tmat2;
        double temp4;
        while(getline(flsr2, line4)){
                vector <double> row;
                istringstream istr1(line4);
                while(istr1 >> temp4){
                        row.push_back(temp4);
                }

                tmat2.push_back(row);
        }
        real_dist = tmat2;*/

      
        real_dist = read_correlation(res, all_size, dist_file, chr_pp);
	simulation.real_dist = real_dist;
	//simulation.online_list = all_list;
        simulation.olist = chr_pp;
	vector <int> x_l;
	vector <int> o_l;
	vector <int> chr_index;
	int olist_count = 0;
	vector <vector <int> > list_start_end;
	for(int i = 0; i < chr_pp.size(); i++){
		vector <int> one_start_end;
		one_start_end.push_back(olist_count);
		//cout<<olist_count;
		olist_count += atoi(chr_pp[i][1].c_str());
		one_start_end.push_back(olist_count);
		list_start_end.push_back(one_start_end);
		//cout<<" "<<olist_count<<endl;
		//for(int j = 0; j < olist2[i].size(); j++){
		//	cout<<olist2[i][j]<<" ";
		//}
		//cout<<endl;
	}
	/*for(int i = 0; i < list_start_end.size(); i++){
		cout<<"list has "<<list_start_end[i][0]<<" "<<list_start_end[i][1]<<endl;
	}*/
	//exit(1);
	//delta2 = calculate_t1(all_size, real_dist, mu2, d0, chr_pp);
        vector <int> loado;
        for(int i = 0; i < chr_pp.size(); i++){
                int i1 = i + 1;
            if(i1 != chr_pp.size()){
                for(int j = 0; j < chr_pp.size(); j++){
                                if(atoi(chr_pp[j][0].c_str()) == i1){
                                        int st = list_start_end[j][0];
                                        int sp = list_start_end[j][1];
                                        for(int j = st; j < sp; j++){
                                                loado.push_back(j);
                                        }
                                }
                }
            }else{
                for(int j = 0; j < chr_pp.size(); j++){
                                if(chr_pp[j][0] == "X"){
                                        int st = list_start_end[j][0];
                                        int sp = list_start_end[j][1];
                                        for(int j = st; j < sp; j++){
                                                loado.push_back(j);
                                        }
                                }
                }
            }
        }
     	char *fls5 = const_cast<char*>(impute_fi.c_str());
        if(impute_fi=="nonei"){
                cout<<"you should insert an imputed hic matrix"<<endl;
                exit(1);
        }
        //char *fls5 = const_cast<char*>(s5.c_str());
        ifstream flsr5;
        flsr5.open(fls5);
        string line5;
        vector <vector<double> > tmat5;
        double temp5;
        while(getline(flsr5, line5)){
                vector <double> row;
                istringstream istr1(line5);
                while(istr1 >> temp5){
                        row.push_back(temp5);
                }

                tmat5.push_back(row);
        }
        for(int i = 0; i < tmat5.size(); i++){
                vector <double> lll;
                for(int j = 0; j < tmat5.size(); j++){
                        lll.push_back(0);
                }
                delta2.push_back(lll);
        }
	/*for(int i = 0; i < loado.size(); i++){
		cout<<i<<" "<<loado[i]<<endl;
	}
	exit(1);*/
        for(int i = 0; i < tmat5.size(); i++){
                for(int j = 0; j < tmat5.size(); j++){
                        int x = loado[i];
                        int y = loado[j];
			//cout<<i<<"	"<<j<<"	should be "<<x<<"	"<<y<<endl;
                        delta2[x][y] = tmat5[i][j];
                }
        }
	simulation.list_start_end = list_start_end;
	vector <vector <vector <int> > > hic_cube;
	vector <vector <int> > plist;
	/*for(int i = 0; i < 655; i++){
		for(int j = 0; j < 655; j++){
			cout<<delta2[i][j]<<" ";
		}
		cout<<endl;
	}*/
	//exit(1);
	/*for(int i = 0; i < olist2.size(); i++){
		for(int j = 0; j < olist2[i][1]; j++){
			vector <int> v1;
			v1.push_back(olist2[i][2]);
			v1.push_back(j);
			plist.push_back(v1);
		}
		
	}*/
	/*for(int i = 0; i < plist.size(); i++){
		cout<<i<<" "<<plist[i][0]<<" "<<plist[i][1]<<endl;
	}*/
	simulation.plist = plist;
	cout<<"ok"<<endl;	
	simulation.hic_cube = hic_cube;
	vector <double> chr_sum1;
	for(int i = 0; i < chr_pp.size(); i++){
		chr_sum1.push_back(0);
	}
	//int r_length = 0;//(biggest_length * 3)/16 + 1;
	simulation.rearrenge_cost = chr_sum1;
	//real_dist2 = real_dist;
	//cout<<"s1"<<endl;
	//real_dist2 = real_dist;
	//delta2 = calculate_t1(all_size, real_dist2, mu2, d0, chr_pp);
        //cout<<"s2"<<endl;
	simulation.delta2 = delta2;
	simulation.online_list.push_back(all_list[0]);
        for(int j = 0; j < simulation.online_list[0].size(); j++){
        	simulation.list2.append_node(simulation.online_list[0][j].x, simulation.online_list[0][j].y, simulation.online_list[0][j].z);
	}        
	for(int k = 1; k <= hic_list.size()-1; k++){
		//simulation.list2.clear();
		cout<<simulation.list2.num_nodes<<endl;
		simulation.online_list.push_back(all_list[k]);
		sim_times = simulation.online_list.size() * 100;
		cout<<sim_times<<endl;
		simulation.real_dist2 = real_dist;
/////////////////////////////////////////////////////////////starting rearrange//////////////////////
        	int failed_times = 0;
		vector <double> num;
		vector <double> energy;
		vector <double> tem;
		temperature = 10;
	       	cout<<"Cooling started."<<endl;

	        for(int i =  1; i <1000000; i++ ){
        	        int counter = 0;
                	int accept_counter = 0;
	
	        	time_t my_time7 = time(NULL);
        		cout<<"@@@@@@@@time@@@@@@@:"<<ctime(&my_time7)<<endl;
                	cout<<"Temperature: "<<temperature<<endl;
                	while(counter < sim_times){
                        	counter++;
		        	timeval t1;
	        		gettimeofday(&t1, NULL);
        			srand(t1.tv_usec * t1.tv_sec);
                	        int rand1 = simulation.online_list.size();//a11_3
                        	if(simulation.simulate_n(k, biggest_length) == true){
                                	accept_counter ++;
                        	}
	                        energy.push_back(simulation.En);
        	                num.push_back(i);
                	        if(accept_counter == sim_times / 10){
                        	        failed_times = 0;
                                	break;
                        	}
                       		if(counter == sim_times - 1){
                                	failed_times++;
                        	}
                	}
			cout<<"Accept number is "<<accept_counter<<endl;
			if(failed_times == 3){
				break;
				//simulation.failed = simulation.failed + 1;
			}
			cout<<"Temperature: "<<temperature<<" finished"<<endl;
			time_t my_time8 = time(NULL);
                	cout<<"@@@@@@@@time@@@@@@@:"<<ctime(&my_time8)<<endl;
	                tem.push_back(temperature);
        	        temperature = temperature * 0.9;
                	simulation.temperature = temperature;
                
        	}

        	cout<<"after rearrange online_list is "<<simulation.online_list.size()<<endl;
                for(int j = 0; j < simulation.online_list[k].size(); j++){
                        simulation.list2.append_node(simulation.online_list[k][j].x, simulation.online_list[k][j].y, simulation.online_list[k][j].z);
                }
		cout<<simulation.list2.num_nodes<<endl;
	}
	ofstream file12;
        file12.open(out_file2);
        for(int i = 0; i < simulation.online_list.size(); i++){
                for(int j = 0; j < simulation.online_list[i].size();j++){
                        if(simulation.online_list[i][j].x2 != 1){
                                file12<<simulation.online_list[i][j].x<<" "<<simulation.online_list[i][j].y<<" "<<simulation.online_list[i][j].z<<endl;;
                        }else{
                                file12<<"!!!!!!0 0 0"<<endl;
                        }
                }
        }
        for(int i = 0; i < chr_pp.size(); i++){
                for(int j = 0; j < chr_pp[i].size(); j++){
                        file12<<chr_pp[i][j]<<" ";
                }
                file12<<endl;
        }
        file12.close();
	        int r_length = 10;
        	//cout<<"get in "<<endl;
	        simulation.no_move(r_length, 10);
	        ofstream file11;
        	file11.open(out_file);
	        for(int i = 0 ; i <simulation.list2.num_nodes; i = i + 1){
        	        file11<<simulation.list2.get_x_value(i+1)/1<<"  "<<simulation.list2.get_y_value(i+1)/1<<"       "<<simulation.list2.get_z_value(i+1)/1<<endl;
	        }
        	file11.close();
        time_t my_time3 = time(NULL);
        cout<<"@@@@@@@@time@@@@@@@:"<<ctime(&my_time3)<<endl;
	return 0;

}

