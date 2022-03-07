#include <iostream.h>
#include <string.h>
#include <stdlib.h>

int rows=-1;
int cols=-1;

double min=0;
double max=1;

int xtics=10;
int xinit=0;
double xinittic=0;
double xfinaltic=0;
int ytics=10;
int yinit=0;
double yinittic=0;
double yfinaltic=0;

int legendtics=5;

double size=6;
double xsize;
double ysize;

char font[20];
int fontsize=16;

char caption_font[20];
int caption_fontsize=16;
char caption_str[1024];

char grey=0;
char revvid=0;
char nolegend=0;

char xlabel[1024];
char xlabel_font[40];
int xlabel_font_size = 16;

char ylabel[1024];
char ylabel_font[40];
int ylabel_font_size = 16;



void car(char* str, char delimiter, char* first)
{
  int i;

  for (i=0; str[i] != delimiter && str[i] != '\0'; i++)
    first[i] = str[i];
  first[i] = '\0';
}

void cdr(char* str, char delimiter, char* rest)
{
  int i,j;

  for (i=0; str[i] != delimiter && str[i] != '\0'; i++);
  for (j=i; str[j] != '\0'; j++)
    rest[j-i] = str[j];
  rest[j-i] = '\0';
}

double trans(double x)
{
  double a = (x - min)/(max-min);
  if (a > 1.0) return 1.0;
  if (a < 0.0) return 0.0;
  return a;
}

double itrans(double y)
{
  return (max-min)*y + min;
}

void psheader(ostream&,double,double,double,double);
void psvariables(ostream&,double,double,int,int,int,char*);
void psdraw(ostream&,istream&,int,int,int,char*,char*,int,char*,char*,int);
void psdrawbw(ostream&,istream&,int,int,int,char*,char*,int,char*,char*,int); 
void pscaption(ostream&,char*,char*,int);
void psclose(ostream&);
void Display_help(ostream&,char*);

main(int argc, char* argv[])
{
  if (argc < 5)
    {
      Display_help(cerr,argv[0]);
      exit (0);
    }

  strcpy(font,"Times-BoldItalic");
  strcpy(caption_font,"Times-BoldItalic");
  strcpy(caption_str," ");
  strcpy(xlabel," ");
  strcpy(ylabel," ");
  strcpy(xlabel_font,"Times-BoldItalic");
  strcpy(ylabel_font,"Times-BoldItalic");

  for (int i=1; i < argc; i++)
    {
      char argstr[50];
      char name[10];
      char tmp[50];
      char val[50];
      strcpy(argstr,argv[i]+2);
      car(argstr,'=',name);
      cdr(argstr,'=',tmp);
      strcpy(val,tmp+1);
      
      if (strcmp(name,"rows") == 0) rows = atoi(val);
      if (strcmp(name,"cols") == 0) cols = atoi(val);

      if (strcmp(name,"min") == 0) min = atof(val);
      if (strcmp(name,"max") == 0) max = atof(val);
      
      if (strcmp(name,"size") == 0) size = atof(val);
      if (strcmp(name,"fontsize") == 0) fontsize = atoi(val);
      if (strcmp(name,"font") == 0) strcpy(font,val);
      if (strcmp(name,"xlabelfontsize") == 0) xlabel_font_size =
						 atoi(val);
      if (strcmp(name,"ylabelfontsize") == 0) ylabel_font_size =
						 atoi(val);
      if (strcmp(name,"captionfontsize") == 0) caption_fontsize =
						 atoi(val);
      if (strcmp(name,"captionfont") == 0) strcpy(caption_font,val);
      if (strcmp(name,"xlabelfont") == 0) strcpy(xlabel_font,val);
      if (strcmp(name,"ylabelfont") == 0) strcpy(ylabel_font,val);
      if (strcmp(name,"caption") == 0) strcpy(caption_str,val);
      if (strcmp(name,"xlabel") == 0) strcpy(xlabel,val);
      if (strcmp(name,"ylabel") == 0) strcpy(ylabel,val);
      if (strcmp(name,"grey") == 0)  grey = 1;
      if (strcmp(name,"gray") == 0)  grey = 1;
      if (strcmp(name,"revvid") == 0)  revvid = 1;
      if (strcmp(name,"nolegend") == 0)  nolegend = 1;

      if (strcmp(name,"xtics") == 0) xtics = atoi(val);
      if (strcmp(name,"ytics") == 0) ytics = atoi(val);
      if (strcmp(name,"xinittic") == 0) xinittic = atof(val);
      if (strcmp(name,"yinittic") == 0) yinittic = atof(val);
      if (strcmp(name,"xfinaltic") == 0) xfinaltic = atof(val);
      if (strcmp(name,"yfinaltic") == 0) yfinaltic = atof(val);


      if (strcmp(name,"legendtics") == 0) legendtics = atoi(val);
      if (strcmp(name,"help") == 0) Display_help(cerr,argv[0]);
    }

  if (rows < 0 || cols < 0) exit(0);

  if (rows == cols)
    xsize = ysize = size;

  if (rows < cols)
    {
      xsize = size;
      ysize = double(rows)/double(cols)*size;
    }

  if (cols < rows)
    {
      ysize = size;
      xsize = double(cols)/double(rows)*size;
    }

  if (rows < ytics) ytics = rows;
  if (cols < xtics) xtics = cols;

  double yminbound= 22.2;
  double ymaxbound = (5.75 + ysize/2) * 72.2+50+caption_fontsize*3;
  double tmpxsize = xsize > 6 ? xsize : 6;
  double xminbound = (4.25 - tmpxsize/2) * 72.2 - 50;
  double xmaxbound = (4.25 + tmpxsize/2) * 72.2 + 50;

  if (xfinaltic == xinittic) xfinaltic = cols;
  if (yfinaltic == yinittic) yfinaltic = rows;

  psheader(cout,xminbound,xmaxbound,yminbound,ymaxbound);
  psvariables(cout,xsize,ysize,rows,cols,fontsize,font);
  if (grey)
    psdrawbw(cout,cin,xtics,ytics,legendtics,xlabel,xlabel_font,xlabel_font_size, ylabel, ylabel_font, ylabel_font_size);
  else
    psdraw(cout,cin,xtics,ytics,legendtics,xlabel,xlabel_font,xlabel_font_size, ylabel, ylabel_font, ylabel_font_size);

  pscaption(cout,caption_str,caption_font,caption_fontsize);
  psclose(cout);
}


void psheader(ostream& output, double xmin, double xmax, double ymin,
	      double ymax)
{
  output << "%!PS-Adobe-2.0\n";
  output << "%%Title: Matrixto.ps\n";
  output << "%%Creator: Sagar A. Pandit \n";
  output << "%%Pages: 1\n";
  output << "%%PageOrder: Ascend\n";
  output << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << 
    " " << ymax << '\n';
  output << "%%DocumentPaperSizes: a4\n";
  output << "%%EndComments\n";
  output << "\n";
  output << "% This is definition for page transformations\n";
  output << "% Do not Change this \n";
  output << "/in {72.2 mul} def\n";
  output << "\n";
  output << "/PageHeight 11.75 in def\n";
  output << "/PageWidth 8.25 in def\n";
  output << "\n";
  output << "/Sat 0.38 def\n";
  output << "%\n";
  output << "% xlen ylen fBox ___\n";
  output << "%\n";
  output << "% where : xlen, ylen are in pt.\n";
  output << "%\n";
  output << "/fBox { \n";
  output << "  /ylen exch def \n";
  output << "  /xlen exch def \n";
  output << "  currentpoint moveto \n";
  output << "  gsave \n";
  output << "    xlen 0 rlineto \n";
  output << "    0 ylen rlineto \n";
  output << "    -1 xlen mul 0 rlineto \n";
  output << "    0 -1 ylen mul rlineto fill \n";
  output << "  grestore \n";
  output << "} def \n";
  output << "/fmod {\n";
  output << " /yy exch def\n";
  output << "  /xx exch def\n";
  output << "  xx xx yy div truncate yy mul sub\n";
  output << "} def\n";
  output << "\n";
  output << "/cmap {\n";
  output << "  /vv exch def\n";
  output << "  1.0 val sub 0.5 add 1.0 fmod val 2.0 mul 1.0 sethsbcolor\n";
  output << "} def\n";
  output << "\n";
  output << "\n";
  output << "%\n";
  output << "% xlen ylen Box ___\n";
  output << "%\n";
  output << "% where : xlen, ylen are in pt.\n";
  output << "%\n";
  output << "/Box { \n";
  output << "  /ylen exch def \n";
  output << "  /xlen exch def \n";
  output << "  currentpoint moveto \n";
  output << "  gsave \n";
  output << "    xlen 0 rlineto \n";
  output << "    0 ylen rlineto \n";
  output << "    -1 xlen mul 0 rlineto \n";
  output << "    0 -1 ylen mul rlineto stroke\n";
  output << "  grestore \n";
  output << "} def \n";
  output << "\n";
  output << "%\n";
  output << "% i j val Cell ___\n";
  output << "%\n";
  output << "% where : i,j are row and column indices starting from 0.\n";
  output << "%              and val is the value of element [0.0:1.0]\n";
  output << "%\n";
  output << "/Cell {\n";
  output << "  /val exch def\n";
  output << "  /j exch def\n";
  output << "  /i exch def\n";
  output << "  \n";
  //  output << "  /val val 2 mul def\n";
  output << "\n";
  output << "  gsave\n";
  output << "    j XCellWidth mul i YCellWidth mul translate\n";
  output << "    0 0 moveto\n";
//   output << "    val 1 gt \n";
//   output << "    { /val val 1 sub def\n";
//   output << "      val 1 val sub 0 setrgbcolor }\n";
//   output << "    { 0 val 1 val sub setrgbcolor } ifelse\n";
//  output << "    val val val 1.0 add sethsbcolor\n";
  output << "    val cmap\n";
  output << "    XCellWidth YCellWidth fBox\n";
  output << "  grestore\n";
  output << "} def\n";
  output << "\n";
  output << "%\n";
  output << "% val BWMap val\n";
  output << "%\n";
  output << "% where : val is the intensity. 0.0 <= val <= 1.0\n";
  output << "%\n\n";
  output << "/BWMap {\n";
  output << " /myval exch def\n";
  output << " myval\n";
  output << "}def\n\n";
  output << "% i j val CellBW ___\n";
  output << "%\n";
  output << "% where : i,j are row and column indices starting from 0.\n";
  output << "%              and val is the value of element [0.0:1.0]\n";
  output << "%\n";
  output << "/CellBW {\n";
  output << "  /val exch def\n";
  output << "  /j exch def\n";
  output << "  /i exch def\n";
  output << "\n";
  output << "  gsave\n";
  output << "    j XCellWidth mul i YCellWidth mul translate\n";
  output << "    0 0 moveto\n";
  if (!revvid)
    output << "    val BWMap setgray\n";
  else
    output << "    1.0 val sub setgray\n";
  output << "    XCellWidth YCellWidth fBox\n";
  output << "  grestore\n";
  output << "} def\n";
  output << "\n";
  output << "\n";
  output << "/XTic {\n";
  output << "  /symbolj exch def\n";
  output << "  /ticj exch def\n";
  output << "  /xticstr 8 string def\n";
  output << "  /xticstr ticj xticstr cvs def\n";
  output << "  gsave\n";
  output << "    0.5 setlinewidth\n";
  output << "    ticj XCellWidth mul 0 translate\n";
  output << "    0 0 moveto\n";
  output << "    0 0 0 setrgbcolor\n";
  output << "    0 -1 TicLen mul rlineto\n";
  output << "    0 -1 5 FontSize add mul rmoveto \n";
  //  output << "    xticstr Cshow stroke\n";
  output << "    symbolj Cshow stroke\n";
  output << "    0 YWidth moveto\n";
  output << "    0 TicLen rlineto stroke\n";
  output << "  grestore\n";
  output << "} def\n";
  output << "\n";
  output << "/YTic {\n";
  output << "  /symboli exch def\n";
  output << "  /tici exch def\n";
  output << "  /yticstr 8 string def\n";
  output << "  /yticstr tici yticstr cvs def\n";
  output << "  gsave\n";
  output << "    0.5 setlinewidth\n";
  output << "    0 tici YCellWidth mul translate\n";
  output << "    0 0 moveto\n";
  output << "    0 0 0 setrgbcolor\n";
  output << "    -1 TicLen mul 0 rlineto \n";
  output << "    -5 0 rmoveto \n";
  //  output << "    yticstr RCshow stroke\n";
  output << "    symboli RCshow stroke\n";
  output << "    XWidth 0 moveto\n";
  output << "    TicLen 0 rlineto stroke\n";
  output << "  grestore\n";
  output << "} def\n";
  output << "\n";
  output << "\n";
  output << "/Legend {\n";
  output << "  0 0 moveto\n";
  output << "  LegendWidth LegendCellWidth add LegendHeight Box\n";
  output << "  gsave\n";
  output << "    0 0 moveto\n";
  output << "    /val 0 def\n";
  output << "    0 1 LegendResolution\n";
  output << "    {\n";
//   output << "      val 1 gt \n";
//   output << "      { val 1 sub 1 val 1 sub sub 0 setrgbcolor }\n";
//   output << "      { 0 val 1 val sub setrgbcolor } ifelse\n";
//  output << "      val val val 1.0 add sethsbcolor\n";
  output << "      val cmap\n";
  output << "      LegendCellWidth LegendHeight fBox\n";
  output << "      /val val 1 LegendResolution div add def\n";
  output << "      LegendCellWidth 0 rmoveto\n";
  output << "    } for\n";
  output << "  grestore\n";
  output << "} def\n";
  output << "\n";
  output << "/LegendBW {\n";
  output << "  0 0 moveto\n";
  output << "  LegendWidth LegendCellWidth add LegendHeight Box\n";
  output << "  gsave\n";
  output << "    0 0 moveto\n";
  output << "    /val 0 def\n";
  output << "    0 1 LegendResolution\n";
  output << "    {\n";
  if (!revvid)
    output << "      val BWMap setgray\n";
  else
    output << "      1.0 val sub setgray\n";
  output << "      LegendCellWidth LegendHeight fBox\n";
  output << "      /val val 1 LegendResolution div add def\n";
  output << "      LegendCellWidth 0 rmoveto\n";
  output << "    } for\n";
  output << "  grestore\n";
  output << "} def\n";
  output << "\n";
  output << "/Cshow { dup stringwidth pop 2 div -1 mul 0 rmoveto show } def\n";
  output << "/RCshow { dup stringwidth -2 div exch  -1 mul exch rmoveto show } def\n";
  output << "\n";
  output << "/LegendTic {\n";
  output << "  /lval exch def\n";
  output << "  /lstr exch def\n";
  output << "  gsave\n";
  output << "    0.5 setlinewidth\n";
  output << "    lval LegendResolution mul LegendCellWidth mul 0 translate\n";
  output << "    0 0 moveto\n";
  output << "    0 0 0 setrgbcolor\n";
  output << "    0 -1 TicLen mul rlineto\n";
  output << "    0 -1 5 FontSize add mul rmoveto\n";
  output << "    lstr Cshow stroke\n";
  output << "  grestore\n";
  output << "} def\n";
  output << "\n";
  output << "%\n";
  output << "% Parameters which change\n";
  output << "%\n";
  output << "\n";
}


void psvariables(ostream& output, double xwidth, double ywidth, int
		 rows, int cols, int fontsize, char* fontname)
{
  output << "/XWidth " << xwidth << " in def\n";
  output << "/YWidth " << ywidth << " in def\n";
  output << "\n";
  output << "/XNumberOfCells " << cols << " def\n";
  output << "/YNumberOfCells " << rows << " def\n";
  output << "\n";
  output << "/TicLen 8 def\n";
  output << "\n";
  output << "/LegendHeight 0.2 in def\n";
  output << "/LegendWidth 6 in def\n";
  output << "/LegendResolution 1000 def\n";
  output << "/LegendYLocation 1 in def\n";
  output << "\n";
  output << "/FontSize " << fontsize << " def\n";
  output << "\n";
  output << "%\n";
  output << "% The global variables \n";
  output << "%\n";
  output << "\n";
  output << "/XLocation PageWidth XWidth sub 2 div def\n";
  output << "/YLocation PageHeight YWidth sub 2 div def\n";
  output << "\n";
  output << "/XCellWidth XWidth XNumberOfCells div def\n";
  output << "/YCellWidth YWidth YNumberOfCells div def\n";
  output << "\n";
  output << "/LegendCellWidth LegendWidth LegendResolution div def\n";
  output << "/LegendXLocation PageWidth LegendWidth sub 2 div def\n";
  output << "\n";
  output << "/" << fontname << " findfont FontSize scalefont setfont\n";
  output << "    \n";
  output << "%\n";
  output << "% Actual drawing\n";
  output << "%\n";
  output << "\n";
}

void psdraw(ostream& output, istream& input, int xtics, int ytics, int
	    legendtics, char* xlabel, char* xlabelfont, int
	    xlabelfontsize, char* ylabel, char* ylabelfont, int
	    ylabelfontsize )
{
  output << "  \n";
  output << "gsave\n";
  output << "  XLocation YLocation translate\n";
  output << "  0 0 moveto\n";
  output << "  0 0 0 setrgbcolor\n";
  output << "  XWidth YWidth Box\n";
  output << "  0.1 setlinewidth\n";
  output << "\n";

 output << " gsave\n";
 output << "  /" << xlabelfont << " findfont " << xlabelfontsize << 
   " scalefont setfont\n"; 
 output << "   XWidth 2 div -1 5 FontSize add " << xlabelfontsize*2 <<
   " add mul moveto \n";
 output << "   (" << xlabel << ") Cshow\n";
 output << " grestore\n";
 output << "\n";
 output << " gsave\n";
 output << "  /" << ylabelfont << " findfont " << ylabelfontsize << 
   " scalefont setfont\n"; 
 output << "   -5 (000) stringwidth pop -1 mul add " <<
   ylabelfontsize*2 << " sub YWidth 2 div moveto \n";
 output << "   90 rotate\n";
 output << "   (" << ylabel << ") Cshow\n";
 output << "   -90 rotate\n";
 output << " grestore\n";

/* for (int i=0; i < rows; i++)
   for (int j=0; j < cols; j++)*/
   int count;
   input >> count;
   for (int i=0; i < count; i++)
     {
       double r,c;
       double v;
       input >> r >> c >> v;
       if ( r <= rows && c <= cols )
         output << "  " << r << " " << c << " " << trans(v) << " Cell ";
       output << '\n';
     }

  for (int i=0; i < xtics+1; i++)
    output << int((cols-xinit)/xtics*i) << " (" <<
      ((xfinaltic-xinittic)/xtics*i)+xinittic << ") " << " XTic\n";
  for (int j=0; j < ytics+1; j++)
    output << int((rows-yinit)/ytics*j) <<  " (" <<
      ((yfinaltic-yinittic)/ytics*j)+yinittic << ") " << " YTic\n";

  int number;
  input >> number;
  cerr << number << '\n';

  for (int i=0; i < number; i++)
    {
      double x1, y1, x2, y2;
      input >> x1 >> y1 >> x2 >> y2; 
      if (x1 <= rows && y1 <= cols) {
      cerr <<  x1 << ' ' <<  y1 << ' ' <<  x2 << ' ' <<  y2 << '\n'; 
      output << "gsave";
      output << "  0 0 0 setrgbcolor\n";
      output << "  1 setlinewidth\n";
      output << "  " << y1 << " XCellWidth mul " << x1 << " YCellWidth mul moveto\n";
      output << "  " << y2 << " XCellWidth mul " << x2 << " YCellWidth mul  lineto\n";
      output << "  stroke\n";
      output << "grestore\n"; }
    }

  output << "grestore\n";
  output << "\n";


  if (!nolegend)
    {
      output << "gsave\n";
      output << "  LegendXLocation LegendYLocation translate\n";
      output << "  Legend\n";
      
      output.width(6);
      for (int k=0; k < legendtics+1; k++)
	output << "(" << itrans(1.0/double(legendtics)*k) << ") " <<
	  1.0/double(legendtics)*k << " LegendTic\n";

      output << "grestore\n";
    }
}

void psdrawbw(ostream& output, istream& input, int xtics, int ytics,
	      int legendtics, char* xlabel, char* xlabelfont, int 
	      xlabelfontsize, char* ylabel, char* ylabelfont, int
	      ylabelfontsize )
{
  output << "  \n";
  output << "gsave\n";
  output << "  XLocation YLocation translate\n";
  output << "  0 0 moveto\n";
  output << "  0 0 0 setrgbcolor\n";
  output << "  XWidth YWidth Box\n";
  output << "  0.1 setlinewidth\n";
  output << "\n";

  output << " gsave\n";
  output << "  /" << xlabelfont << " findfont " << xlabelfontsize << 
    " scalefont setfont\n"; 
  output << "   XWidth 2 div -1 5 FontSize add " << xlabelfontsize*2 <<
    " add mul moveto \n";
  output << "   (" << xlabel << ") Cshow\n";
  output << " grestore\n";
  output << "\n";
  output << " gsave\n";
 output << "  /" << ylabelfont << " findfont " << ylabelfontsize << 
   " scalefont setfont\n"; 
 output << "   -5 (000) stringwidth pop -1 mul add " <<
   ylabelfontsize*2 << " sub YWidth 2 div moveto \n";
 output << "   90 rotate\n";
 output << "   (" << ylabel << ") Cshow\n";
 output << "   -90 rotate\n";
 output << " grestore\n";


/* for (int i=0; i < rows; i++)
   for (int j=0; j < cols; j++) */
   int count;
   input >> count;
   for (int i=0; i < count; i++)
       { 
	 double r,c;
	 double v;
	 input >> r >> c >> v;
         if ( r <= rows && c <= cols )
	   output << "  " << r << " " << c << " " << trans(v) << " CellBW ";
	 output << '\n';
       }

  for (int i=0; i < xtics+1; i++)
    output << int((cols-xinit)/xtics*i) << " (" <<
      ((xfinaltic-xinittic)/xtics*i) + xinittic << ") " << " XTic\n";
  for (int j=0; j < ytics+1; j++)
    output << int((rows-yinit)/ytics*j) <<  " (" <<
      ((yfinaltic-yinittic)/ytics*j) + yinittic << ") " << " YTic\n";

  int number;
  input >> number;
  cerr << number << '\n';

  for (int i=0; i < number; i++)
    {
      double x1, y1, x2, y2;
      input >> x1 >> y1 >> x2 >> y2; 
      if (x1 <= rows && y1 <= cols) {
      cerr <<  x1 << ' ' <<  y1 << ' ' <<  x2 << ' ' <<  y2 << '\n'; 
      output << "gsave";
      output << "  0 0 0 setrgbcolor\n";
      output << "  3 setlinewidth\n";
      output << "  " << y1 << " XCellWidth mul " << x1 << " YCellWidth mul moveto\n";
      output << "  " << y2 << " XCellWidth mul " << x2 << " YCellWidth mul  lineto\n";
      output << "  stroke\n";
      output << "grestore\n";}
    }

  output << "grestore\n";
  output << "\n";
  if (!nolegend)
    {
      output << "gsave\n";
      output << "  LegendXLocation LegendYLocation translate\n";
      output << "  LegendBW\n";
      
      output.width(6);
      for (int k=0; k < legendtics+1; k++)
	output << "(" << itrans(1.0/double(legendtics)*k) << ") " <<
	  1.0/double(legendtics)*k << " LegendTic\n";

      output << "grestore\n";
    }
}

void pscaption(ostream& output, char* caption_str, char* caption_font, 
	       int caption_fontsize)
{
  output << "clear\n";
  output << "/" << caption_font << " findfont " << caption_fontsize << 
    " scalefont setfont\n";
  output << "PageWidth 2 div YLocation YWidth add " <<
    1.8*caption_fontsize << " add moveto\n"; 
  output << "(" << caption_str << ") Cshow\n";
  output << "stroke\n";
}

void psclose(ostream& output)
{
  output << "showpage\n";
}

void Display_help(ostream& error, char* prog_name)
{
  error << "Usage: " << prog_name << " <Mandatory options> <Other options>\n";
  error << "\n";
  error << "Reads the matrix on STDIN and writes Postscript file on STDOUT.\n";
  error << "Format of the input matrix is: <i> <j> <value>\n";
  //error << "\n";
  error << "Mandatory options:\n";
  error << "  --rows=<Number of rows in the matrix>          [No default value]\n";
  error << "  --cols=<Number of columns in the matrix>       [No default value]\n";
  error << "  --min=<Minimum element in the matrix>          [Default 0.0]\n";
  error << "  --max=<Maximum element in the matrix>          [Default 1.0]\n";
  //error << "\n";
  error << "Other options:\n";
  error << "  --xtics=<Number of tics on column index>       [Default 10]\n";
  error << "  --ytics=<Number of tics on row index>          [Default 10]\n";

  error << "  --xinittic=<Initial tic on column index>       [Default 0]\n";
  error << "  --yinittic=<Initial tic on row index>          [Default 0]\n";
  error << "  --xfinaltic=<Final tic on column index>        [Default 1]\n";
  error << "  --yfinaltic=<Final tic on row index>           [Default 1]\n";

  error << "  --legendtics=<Number of tics on legend bar>    [Default 5]\n";
  error << "  --font=<Valid Postscript font name>            [Default Times-BoldItalic]\n";
  error << "  --fontsize=<Size of the font>                  [Default 16]\n";
  error << "  --caption=<String quoted in \"\">                [Default \"\"]\n";
  error << "  --captionfont=<Valid Postscript font name>     [Default Times-BoldItalic]\n";
  error << "  --captionfontsize=<Size of the font>           [Default 16]\n";
  error << "  --xlabel=<String quoted in \"\">               [Default \"\"]\n";
  error << "  --xlabelfont=<Valid Postscript font name>      [Default Times-BoldItalic]\n";
  error << "  --xlabelfontsize=<Size of the font>            [Default 16]\n";
  error << "  --ylabel=<String quoted in \"\">               [Default \"\"]\n";
  error << "  --ylabelfont=<Valid Postscript font name>      [Default Times-BoldItalic]\n";
  error << "  --ylabelfontsize=<Size of the font>            [Default 16]\n";
  error << "  --grey\n";
  error << "  --gray\n";
  error << "  --revvid\n";
  error << "  --nolegend\n";
  //error << "\n";
  error << "e.g.: " << prog_name << " --rows=100 --cols=100 --min=-5 --max=23.6 --gray < The_Matrix\n";
  error << "Converts a 100x100 matrix with values between -5:23.6 to a gray scale image\n";
  error << "Note: There is no space between = and the values.\n";
  exit(0);
} 
