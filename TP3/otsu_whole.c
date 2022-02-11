#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
int **data;
int hist[255];
double sigma[255],sort[255];
void writePGM(int height,int width,int maxVal);
void getPGM(char filename[]);
int histogram_otsu(int,int);
double otsu(int,int,int,int);
void binary_conversion(int t,int ht,int wth,int maxVal);
/**********Function to read a ppm image file*********/
void getPGM(char filename[])
{

	//Assignment operations
	FILE *ifp;	
	int i,j,type,width,height,maxVal,t;	
	int ch_int,ch,n;
	//char ch,r,g,b,hi,lo,R,G,B;
	char file_name[11];
	printf("\n*******************************\n");
	printf("\nOpening the file %s\n",filename);
	printf("\n*******************************\n");
	n= strlen(filename);
	
	i=0;
	while(1)
	{
		if(filename[i]=='\n')
		{
			filename[i]='\0';
			break;
		}
		i++;
	}
	
	//printf("Reading the image file:%s\n",filename); 
	 
	 printf("\n The number of character :%d\n",n);
		
	//File operation for input file
	ifp=fopen(filename,"rb");
	if(ifp==NULL)
	{
		printf("Error: Unable to open file %s\n\n", filename);
		exit(1);
	}
	printf("Reading the image file:%s\n",filename);
	
	
	//Reading data	
	ch=getc(ifp);
	if(ch!='P')
	{
		printf("Not a ppm file...\n");
		exit(1);
	}
	
	ch=getc(ifp);
	putchar(ch);
	type=ch-48;
	
	
	if(type==2)
	{
		printf("It is a P2 file...\n");
	}
	if(type==5)
	{
		printf("It is a P5 file...\n");
	}
	
	
	while((ch=getc(ifp))!='\n');
	ch=getc(ifp);
	

	if(ch=='#')
	{
		
		while(getc(ifp)!='\n');//skiping the comment
		//exit(1);
	}
	fseek(ifp,-1,SEEK_CUR);
	
	
	//Reading width and writing width to output file
	fscanf(ifp,"%d",&width);	
	
	
	//Reading height and writing width to output file	
	fscanf(ifp,"%d",&height);
	
	
	
	//Decalring the dynamic array of type integer for RED, GREEN and BLUE component
	data = (int **)malloc(width*sizeof(int *));
	for(i=0;i<height;i++)
		data[i]=(int *)malloc(height*sizeof(int *));
	
	
	//Extracting the maximum gray value and writing to output file
	fscanf(ifp,"%d\n",&maxVal);
	printf("%d\n",maxVal);
	
	
	
	
	
	if(type==5)//for P6 image file
	{
		
        	for (i = 0; i < height; ++i)
        	{
        		for (j = 0; j < width; ++j) 
        		{

					fread((void *)&data[i][j],sizeof(unsigned char),1,ifp);
					//printf("%d ",data[i][j]);        		
            }
        }
       
      }
          
      /*for(i=0;i<height;i++)
      {
      	for(j=0;j<width;j++)
      	{
      			printf("%d ",data[i][j]);
      	}
      }*/
      
      
      t = histogram_otsu(height,width); 
      binary_conversion(t,height,width,255);  	
     printf("\nThe Threshold ---------------%d",t);
      
	fclose(ifp);
	
}



int histogram_otsu(int ht,int wth)
{

	int i,j,k,t,count,final_threshold;
	float threshold;
	double wt=0.0,mean=0.0,wt_bg,sum=0.0,mean_bg,var=0.0,var_bg;
	double wt1=0.0,mean1=0.0,wt_fg,sum1=0.0,mean_fg,var1=0.0,var_fg;
	double thes,temp;
	
	for(k=0;k<=255;k++)
	{
		count = 0;
		
		for(i=0;i<ht;i++)
		{
			for(j=0;j<wth;j++)
			{
				if(data[i][j]==k)
				{
					count++;
					hist[k]=count;
				}
			}
		}
	}
	
	for(i=0;i<=255;i++)
	{
		printf("%d level total gray values = %d \n",i,hist[i]);
	}
	
	
	//calling OTSU
	for(t=0;t<254;t++)
	{
		//printf("\nLower boundary 0, higher boundary %d",t);
		var_bg = otsu(0,t,ht,wth);
		//printf("\nLower boundary %d, higher boundary 255",t);
		var_fg = otsu(t+1, 255,ht,wth);
		
		//printf("\n The variance background =%f varience forgraound=%f",var_bg,var_fg);
		
		sigma[t] = var_bg + var_fg;	
		sort[t]= sigma[t];	
		
	}
	
	
	//printf("\nThe minimum value is %f",thes);
	for (i=1;i<254;i++)
		{
		//printf("%f ",sigma[i]);
			for(j=1;j<254-i;j++)
			{
					
				if(sigma[i]<sigma[j])
				{
					temp = sigma[i];
					sigma[i]=sigma[j];
					sigma[j]=temp;
				}
			}
		}
		printf("\n");
		/*for(i=1;i<254;i++)
		{
			printf("%f ",sigma[i]);
		}*/			
	thes = sigma[2];
	//printf("Threshold :%f",thes);
	for(i=1;i<=255;i++)
	{
		if(sort[i] == thes)
			threshold = i;
	}
	final_threshold = (int)threshold;
	//printf("\nThe Threshold is %d\n",final_threshold);
	return (final_threshold);
	
	//binary_conversion(final_threshold,ht,wth,255);

}

double otsu(int min,int max,int ht,int wth)
{
	//printf("\n Function call ");
	int i,j,k,t,count,threshold;
	double wt=0.0,mean=0.0,wt_bg,sum=0.0,mean_bg,var=0.0,var_bg;
	double wt1=0.0,mean1=0.0,wt_fg,sum1=0.0,mean_fg,var1=0.0,var_fg;
	double thes,temp,varience,sq,diff;
	
		for(i=min;i<=max;i++)
		{
			//wt = wt + hist[i];
			mean=mean+i*hist[i];
			sum = sum+hist[i];
		}
		if(sum!=0.0)
		{
		//printf("\nMean=%f Sum=%f",mean,sum);
		wt_bg = sum/(ht*wth);
		//printf("\nWeight = %f",wt_bg);
		mean_bg = mean/sum;		
		//printf("\nmean = %f",mean_bg);
		for(j=min;j<=max;j++)
		{
			diff = j-mean_bg;
			sq= pow(diff,2);
			var =var + (sq*hist[j]);
		}
		//printf("\n%f", var);
		var_bg=var/sum;
		varience = wt_bg * pow(var_bg,2);
		return (varience);
		}
		else
			return(0.0);
}
void binary_conversion(int t,int ht,int wth,int maxval)
{
	FILE *output;
	int i,j;
	char ot[11];
	static int c;
	
	sprintf(ot, "otsu%d.pgm",c);c++;
	
	
	
	output = fopen(ot,"wb");
	
	fprintf(output,"P5\n");
	fprintf(output,"%d ",wth);
	fprintf(output,"%d\n",ht);
	fprintf(output,"%d\n",maxval);
	for(i=0;i<ht;i++)
	{
		for(j=0;j<wth;j++)
		{
		
			//printf("%d ",data[i][j]);
			if(data[i][j]<t)
			{
				data[i][j]=0;
			}
			else
			{
				data[i][j]=255;
			}
		}
	}
	
	for(i=0;i<ht;i++)
	{
		for(j=0;j<wth;j++)
		{
			
			fwrite((void *)&data[i][j],sizeof(unsigned char),1,output);
			
		}
	}
	fclose(output);	
}


/************MAIN Function**************************/
main()
{
	int i,j;
	FILE *ip;
	
	ip=fopen("cyan_image_list.txt","rb");
  	char buf[50];
  	while (fgets (buf, sizeof(buf), ip)) 
  	{
    	printf("\nTrying to read file : %s", buf);
   		 //store(buf);
   		getPGM(buf);
  	}
  
	if (ferror(stdin)) 
	{
    fprintf(stderr,"Oops, error reading stdin\n");
    abort();
  }
  return 0;
}
