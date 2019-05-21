//compile: gcc -o parser parser.c
//ref: http://stackoverflow.com/questions/3501338/c-read-file-line-by-line
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset */
#include <unistd.h> /* close */

//src and dst are in two string
struct src_dst_rtt{
	char s_ip[16];
	char d_ip[16];
	//char rtt[10];
	long rtt;
};

//src and dst are in one string
struct srcdst_rtt{
	char sd_ip[32]; //combine src and dst together
	//char rtt[10];
	long rtt;
};

//for struct srcdst_rtt
void parse_1(FILE *from, FILE *to){
	struct src_dst_rtt sdr;
	char s_ip[16];
	char d_ip[16];
	char rtt[10];
	
	char * line = NULL; //each line
	size_t len = 0; //len of each line
	int n; // position in each line
	int i; // # of ','
	char * p = NULL; //position pointer in each line
	ssize_t read; //return size of each line
	while ((read = getline(&line, &len, from)) != -1) {
		//printf("Retrieved line of length %zu :\n", read);
		//printf("%s", line);
	  	//printf("len is %d\n", len);
		p = line;
		i = 0; //position of ","
		n = 0; //# of char
		//printf("line size is: %d\n", strlen(line));
		
		memset(s_ip, 0, 16);
		memset(d_ip, 0, 16);
		memset(rtt, 0, 10);
		
		while(n < read){
			if(i == 3){
				int j = 0;
				while(*p != ','){
					s_ip[j] = *p;
					j++;
					p++;
					n++;
				}
				//printf("add %d to s_ip\n",j);
			}else if(i == 6){
				int j = 0;
				while(*p != ','){
					d_ip[j] = *p;
					j++;
					p++;
					n++;
				}
				//printf("add %d to s_ip\n",j);
			}else if(i == 9){
				int j = 0;
				while(*p != ','){
					rtt[j] = *p;
					j++;
					p++;
					n++;
				}
				//printf("add %d to s_ip\n",j);
			}
			if(*p == ','){
				i++;
				p++;
				n++;
			}else{
				p++;
				n++;
			}
		}
		//printf("s_ip is %s, d_ip is %s, rtt is %s\n", s_ip, d_ip, rtt);
		//ref: http://stackoverflow.com/questions/16645583/how-to-copy-char-array-to-another-char-array-in-c
		strncpy(sdr.s_ip , s_ip, 16); 
		strncpy(sdr.d_ip , d_ip, 16); 
		//strncpy(sdr.rtt , rtt, 10); 
		sdr.rtt = atoi(rtt);
		//printf("%s, %s, %ld\n", sdr.s_ip, sdr.d_ip, sdr.rtt);
		//printf("i = atoi (num) is %d\n", atoi(sdr.rtt));
		//ref: http://stackoverflow.com/questions/15643870/read-write-structures-to-file-c
		//fwrite(Data[i], sizeof (Student), 1, file);
		if(*s_ip != '\0'){
			fwrite(&sdr, sizeof(struct src_dst_rtt), 1, to);
		}
	}
}

// for struct srcdst_rtt
void read_parse_result_1(FILE *from_file){
	struct src_dst_rtt sdr;
	//ref: http://stackoverflow.com/questions/15643870/read-write-structures-to-file-c
	//fread(&Data, sizeof(Student), count, file);
	while(fread(&sdr, sizeof(struct src_dst_rtt), 1, from_file)){
		printf("%s, %s, %ld\n", sdr.s_ip, sdr.d_ip, sdr.rtt);
		memset(&sdr, 0, sizeof(struct src_dst_rtt));
	}
}

//for struct src_dst_rtt
void parse_2(FILE *from, FILE *to){
	struct srcdst_rtt sdr;
	//char s_ip[16];
	//char d_ip[16];
	char sd_ip[32];
	char rtt[10];
	
	char * line = NULL; //each line
	size_t len = 0; //len of each line
	int n; // position in each line
	int i; // # of ','
	char * p = NULL; //position pointer in each line
	ssize_t read; //return size of each line
	int j;
	while ((read = getline(&line, &len, from)) != -1) {
		//printf("Retrieved line of length %zu :\n", read);
		//printf("%s", line);
	  	//printf("len is %d\n", len);
		p = line;
		i = 0; //position of ","
		n = 0; //# of char
		//printf("line size is: %d\n", strlen(line));
		
		memset(sd_ip, 0, 32);
		memset(rtt, 0, 10);
		
		while(n < read){
			if(i == 3){
				j = 0;
				while(*p != ','){
					sd_ip[j] = *p;
					j++;
					p++;
					n++;
				}
				//printf("add %d to s_ip\n",j);
				sd_ip[j] = '-'; //10.28.76.46-10.23.45.67
				j++;
			}else if(i == 6){
				//int j = 0;
				while(*p != ','){
					sd_ip[j] = *p;
					j++;
					p++;
					n++;
				}
				//printf("add %d to s_ip\n",j);
			}else if(i == 9){
				j = 0;
				while(*p != ','){
					rtt[j] = *p;
					j++;
					p++;
					n++;
				}
				//printf("add %d to s_ip\n",j);
			}
			if(*p == ','){
				i++;
				p++;
				n++;
			}else{
				p++;
				n++;
			}
		}
		//printf("s_ip is %s, d_ip is %s, rtt is %s\n", s_ip, d_ip, rtt);
		//ref: http://stackoverflow.com/questions/16645583/how-to-copy-char-array-to-another-char-array-in-c
		strncpy(sdr.sd_ip , sd_ip, 32); 
		//strncpy(sdr.d_ip , d_ip, 16); 
		//strncpy(sdr.rtt , rtt, 10); 
		sdr.rtt = atoi(rtt);
		//printf("%s, %s, %ld\n", sdr.s_ip, sdr.d_ip, sdr.rtt);
		//printf("i = atoi (num) is %d\n", atoi(sdr.rtt));
		//ref: http://stackoverflow.com/questions/15643870/read-write-structures-to-file-c
		//fwrite(Data[i], sizeof (Student), 1, file);
		if(*sd_ip != '\0'){
			fwrite(&sdr, sizeof(struct src_dst_rtt), 1, to);
		}
	}
}

// for struct src_dst_rtt
void read_parse_result_2(FILE *from_file){
	struct srcdst_rtt sdr;
	//ref: http://stackoverflow.com/questions/15643870/read-write-structures-to-file-c
	//fread(&Data, sizeof(Student), count, file);
	long num = 0;
	while(fread(&sdr, sizeof(struct src_dst_rtt), 1, from_file)){
		printf("%s, %ld\n", sdr.sd_ip, sdr.rtt);
		memset(&sdr, 0, sizeof(struct src_dst_rtt));
		num ++;
	}
	printf("%ld records in total\n", num);
}

// test code
void get_file_name(){
	char * str0 = "/tmp/region_Raw_PingmeshData_0";
	char * str1 = "/tmp/region_Raw_PingmeshData_1";
	char fnames[80];
	char num[2];
	int i;
	for(i = 22; i < 100; i++){
		//itoa(i, num, 2);
		sprintf(num, "%d", i);
		//printf("%s \n", num);
		memset(fnames, 0, 80);		
		strcat(fnames, str0);
		printf("%s\n", strcat(strcat(fnames, num),".log"));
	}
}

//parse all input files, and write results to one file
void collect_data(){
	FILE * fp_from;
	FILE * fp_to;
	char * str0 = "/shared/ATC17/PingmeshData/region_Raw_PingmeshData_0";
	char * str1 = "/shared/ATC17/PingmeshData/region_Raw_PingmeshData_";
	char fnames[80];
	char num[2];
	int i;
	char num3[3];

	fp_to = fopen("/tmp/region_Raw_PingmeshData.result", "ab");
	if (fp_to == NULL){
		exit(EXIT_FAILURE);
	}

	for(i = 22; i < 100; i++){
	//for(i = 22; i < 24; i++){
		//itoa(i, num, 2);
		sprintf(num, "%d", i);
		//printf("%s \n", num);
		memset(fnames, 0, 80);		
		strcat(fnames, str0);
		//printf("%s\n", strcat(strcat(fnames, num),".log"));

		fp_from = fopen(strcat(strcat(fnames, num),".log"), "r");
		if(fp_from == NULL){
			exit(EXIT_FAILURE);
		}
		printf("parsing: %s\n", fnames);
		parse_2(fp_from, fp_to);
		fclose(fp_from);
	}
	
	for(i = 100; i < 122; i++){
	//for(i = 22; i < 24; i++){
		//itoa(i, num, 2);
		sprintf(num3, "%d", i);
		printf("%s \n", num3);
		memset(fnames, 0, 80);		
		strcat(fnames, str1);
		//printf("%s\n", strcat(strcat(fnames, num),".log"));

		fp_from = fopen(strcat(strcat(fnames, num3),".log"), "r");
		if(fp_from == NULL){
			exit(EXIT_FAILURE);
		}
		printf("parsing: %s\n", fnames);
		parse_2(fp_from, fp_to);
		fclose(fp_from);
	}
	
	
	fclose(fp_to);
	
	fp_to = fopen("/tmp/region_Raw_PingmeshData.result", "rb");
	if (fp_to == NULL){
		exit(EXIT_FAILURE);
	}
	printf("-----------------------------------------\n");
	read_parse_result_2(fp_to);
	fclose(fp_to);	

}



int main(void)
{
	collect_data();
	//get_file_name();
/*
	FILE * fp_from;
	FILE * fp_to;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	char *p = NULL;
	int i;
	int n = 0;

	char s_ip[16];
	char d_ip[16];
	char rtt[10];

	//start: test 1
	//struct src_dst_rtt sdr;
	fp_from = fopen("/tmp/region_Raw_PingmeshData_022.log", "r");
	if (fp_from == NULL){
		exit(EXIT_FAILURE);
	}
	fp_to = fopen("/tmp/region_Raw_PingmeshData_022.parse1", "wb");
	if (fp_to == NULL){
		exit(EXIT_FAILURE);
	}
	parse_1(fp_from, fp_to);
	fclose(fp_from);
	fclose(fp_to);
	if (line){
		free(line);
	}

	fp_to = fopen("/tmp/region_Raw_PingmeshData_022.parse1", "rb");
	if (fp_to == NULL){
		exit(EXIT_FAILURE);
	}
	printf("-----------------------------------------\n");
	read_parse_result_1(fp_to);
	fclose(fp_to);	
	//end test 1	
	
	//start: test 2
	//struct srcdst_rtt sdr;
	fp_from = fopen("/tmp/region_Raw_PingmeshData_022.log", "r");
	if (fp_from == NULL){
		exit(EXIT_FAILURE);
	}
	fp_to = fopen("/tmp/region_Raw_PingmeshData_022.parse2", "wb");
	if (fp_to == NULL){
		exit(EXIT_FAILURE);
	}
	parse_2(fp_from, fp_to);
	fclose(fp_from);
	fclose(fp_to);
	if (line){
		free(line);
	}

	fp_to = fopen("/tmp/region_Raw_PingmeshData_022.parse2", "rb");
	if (fp_to == NULL){
		exit(EXIT_FAILURE);
	}
	printf("-----------------------------------------\n");
	read_parse_result_2(fp_to);
	fclose(fp_to);	
	//end test 2	
*/	
	exit(EXIT_SUCCESS);
}
