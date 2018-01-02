#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include<time.h>

using namespace std;



vector<string > guess_board;



bool eliminate(int row_id,int col_id,char value,vector< vector< string > > &row_pattern, vector< vector<string > > &col_pattern)
{
    bool change=false;
    vector<string>::iterator it=row_pattern[row_id].begin();
    while(it!=row_pattern[row_id].end())
    {
        if((*it)[col_id]!=value)
        {
       
            row_pattern[row_id].erase(it);
            
        }
        else
        {
            it++;
        }
    }
    it=col_pattern[col_id].begin();
    while(it!=col_pattern[col_id].end())
    {
        if((*it)[row_id]!=value)
        {
    
            col_pattern[col_id].erase(it);
            change=true;
        }
        else
        {
            it++;
        }
    }
    return change;
}

bool check(int &n,vector<vector<int> > &col_hint,vector<string > board)//check if the board dismatch the col-hint
{
    int num=board.size();
    int leftspace=n-num;
    for(int i=0;i<n;i++)
    {
        int ptr=0;
        int count=0;
        int idx_hint=0;
        while(ptr<num)
        {
            if(board[ptr][i]=='1')
            {
                count++;
                ptr++;
                if(idx_hint>=col_hint[i].size()){
                    // cout<<"danger"<<endl;
                    return false;
                }
                else if(count>col_hint[i][idx_hint])
                {
                    // cout<<"hereeliminate"<<endl;
                    return false;
                }
                if(ptr==num)
                {
                    
                    if ((num==n)&&(count!=col_hint[i][idx_hint]))
                    {
                        // cout<<"danger"<<endl;
                        return false;
                    }
                    
                    
                }
            }
            else
            {
                if(count!=0)
                {
                    if(idx_hint>=col_hint[i].size()){
                        // cout<<"danger"<<endl;
                        return false;
                    }
                    if(count==col_hint[i][idx_hint])
                    {
                        idx_hint++;
                    }
                    else
                    {
                        // cout<<"danger"<<endl;
                        return false;
                    }
                    count=0;
                }
                ptr++;
            }
            if(ptr==num)
            {
                int needspace=0;
                int betweenspace=0;
                for(int j=idx_hint;j<col_hint[i].size();j++)
                {
                    needspace+=col_hint[i][j];
                    betweenspace++;
                }
                if(count+leftspace<needspace+betweenspace-1)
                {
                    // cout<<"eliminate in check"<<endl;
                    return false;
                }
            }
            
        }
    }
   
    if(num<n){
        return true;
    }
    // cout<<"success!!"<<endl;
    // for(int i=0;i<n;i++)
    // {
    //     cout<<board[i]<<endl;
    // }
    return true;
}

//idx_row後依據col的生成board
void generateguess(bool fromcol,int idx_row,int n,vector<string >&board,vector <vector<string > > &row_pattern,vector <vector<string > > &col_pattern)
{
    if(fromcol)
    {
        for(int i=idx_row+1;i<n;i++)
        {
        size_t j =board[i].find_first_of('*',0);
        while(j!=string::npos)
        {
            bool change=true;
            vector<string >::iterator it=col_pattern[j].begin();
            char tmp=(*it)[i];
            while(it!=col_pattern[j].end())
            {
                
                if(tmp!=(*it)[i])
                {
                    change=false;
                    break;
                }
                it++;
            }
            if(change)
            {
                board[i][j]=tmp;
            }
            j=board[i].find_first_of('*',j+1);
        }
        
        }
    }
    else
    {
        for(int i=idx_row+1;i<n;i++)
        {
        size_t j =board[i].find_first_of('*',0);
        while(j!=string::npos)
        {
            bool change=true;
            vector<string >::iterator it=row_pattern[i].begin();
            char tmp=(*it)[j];
            while(it!=row_pattern[i].end())
            {
                
                if(tmp!=(*it)[j])
                {
                    change=false;
                    break;
                }
                it++;
            }
            if(change)
            {
                board[i][j]=tmp;
            }
            j=board[i].find_first_of('*',j+1);
        }
        
        }

    }

}


bool preprocess(int idx_row,int n,vector<string > &boardtmp,vector <vector<string > > &row_pattern,vector <vector<string > > &col_pattern)
{
    // cout<<"prepos"<<endl;
    if(idx_row<0)
    {
        return true;
    }

    vector<string > board;
    for(int i=0;i<boardtmp.size();i++)
    {
        board.push_back(boardtmp[i]);
    }
    
    string num;
    for(int i=0;i<n;i++)
    {
        num.push_back('*');
    }
    for(int i=idx_row+1;i<n;i++)
    {
        board.push_back(num);
    }

    //for the idx_row(painted), remove the unmatched col_pattern
    for(int i=0;i<n;i++)
    {
        vector<string >::iterator it=col_pattern[i].begin();
        // cout<<"board:"<<idx_row<<","<<i<<endl;
        while(it!=col_pattern[i].end())
        {
            if((*it)[idx_row]!=board[idx_row][i])
            {
                col_pattern[i].erase(it);
                if(col_pattern[i].empty())
                {
                    // cout<<"colempty"<<i<<endl;
                    return false;
                }
            }
            else
            {
                it++;
            }
        } 
    }

    bool change=true;
    // cout<<"beforeguess"<<endl;
    while(change)
    {
    change=false;
    generateguess(true,idx_row,n,board,row_pattern,col_pattern);

    // cout<<"afterguess"<<endl;
    // 依據board消除row_pattern
    for(int i=idx_row+1;i<n;i++)
    {
        
        size_t j =board[i].find_first_not_of('*',0);
        while(j!=string::npos)
        {
            vector<string >::iterator it=row_pattern[i].begin();
            while(it!=row_pattern[i].end())
            {

                // cout<<"erasewrong"<<endl;   
                if((*it)[j]!=board[i][j])
                {
                    row_pattern[i].erase(it);
                    change=true;
                    if(row_pattern[i].empty())
                    {
                        return false;
                    }
                    continue;
                }
                it++;
            }
            j=board[i].find_first_not_of('*',j+1);
        }
    }
    generateguess(false,idx_row,n,board,row_pattern,col_pattern);

    for(int i=idx_row+1;i<n;i++)
    {
        
        size_t j =board[i].find_first_not_of('*',0);
        while(j!=string::npos)
        {
            vector<string >::iterator it=col_pattern[j].begin();
            while(it!=col_pattern[j].end())
            {

                // cout<<"erasewrong"<<endl;   
                if((*it)[i]!=board[i][j])
                {
                    col_pattern[j].erase(it);
                    change=true;
                    if(col_pattern[j].empty())
                    {
                        return false;
                    }
                    continue;
                }
                it++;
            }
            j=board[i].find_first_not_of('*',j+1);
        }

    }
    }
    
    return true;
    
}


bool DFS(int idx_row,int n,vector<vector<int> > &col_hint,vector<string > board,vector <vector<string > > &row_pattern,vector <vector<string > > &col_pattern,string &result)
{
    vector<vector<string > > new_row_pattern;
    vector<vector<string > > new_col_pattern;

    // cout<<"DFS,idx_row:"<<idx_row<<endl;
    for(int i=0;i<n;i++)
    {
        new_row_pattern.push_back(row_pattern[i]);
        new_col_pattern.push_back(col_pattern[i]);
    }
    
    // cout<<"idx_row:"<<idx_row<<endl;
    if(!check(n,col_hint,board))
    {
        // cout<<"check=true"<<endl;
        return false;
    }
    else if (board.size()==n)//
    {
        cout<<"successfile!!"<<endl;
        for(int i=0;i<n;i++)
        {
            result+=board[i][0];
            for(int j=1;j<n;j++)
            {
                result+="\t";
                result+=board[i][j];
                
            }
            result+="\n";
        }
        return true;
    }
    
    if(!preprocess(idx_row-1,n,board,new_row_pattern,new_col_pattern))
    {
        
        // cout<<"false:"<<endl;
        return false;
    }
   

    
    bool ans;
    for(int i=0;i<new_row_pattern[idx_row].size();i++)//i will iterate the possible pattern
    {
        // cout<<"pattern:"<<i<<endl;
        // cout<<new_row_pattern[idx_row][i]<<endl;
        // for(int j=0;j<n;j++)
        // {
        //     cout<<new_row_pattern[idx_row][i][j];
        // }
        // cout<<endl;
        board.push_back(new_row_pattern[idx_row][i]);
        ans=DFS(idx_row+1,n,col_hint,board, new_row_pattern,new_col_pattern,result);
        if(ans==true)
        {
            break;
        }
        board.pop_back();
    }
    // if(ans==false)
    // {
    //     cout<<"falsefalse"<<endl;
    // }
    return ans;
}


void generate_row(int idx_row,int current,int n,vector<int>&row_hint, string pattern, vector< vector<string > > &row_pattern)//according to the row hint to generate_row
{
    
    if (row_hint.empty())
    {
        // cout<<"empty"<<endl;
        string tmp="";
        for(int i =0;i<n;i++)
        {
            tmp+='0';
        }
        row_pattern[idx_row].push_back(tmp);
        for(int i=0;i<n;i++)
        {
            guess_board[idx_row][i]='0';
        }
        return;
    }
    else if(current>=row_hint.size())
    {
        int space=n-pattern.size();
        if(space<0)
        {
            pattern.pop_back();
        }
        else
        {
            for(int i=0;i<space;i++)
            {
                pattern.push_back('0');
            }
        }
        
        row_pattern[idx_row].push_back(pattern);
        return;
    }
    else
    {   
        // cout<<"current:"<<current<<endl;
        int sum=0;//sum: # need to be 1
        for(int i=current;i<row_hint.size();i++)
        {
            sum+=row_hint[i];
        }
        int space=n-pattern.size()-sum+row_hint.size()-current-1;//row_hint-current-1 : necessary space between left pattern 
        // cout<<"space:"<<space<<endl;
        if(space<0){return;}
        for(int i=0;i<=space;i++)
        {
            string tmp=pattern;
            
            for(int j=0;j<row_hint[current];j++)
            {
                tmp.push_back('1');
            }
            tmp.push_back('0');
            generate_row(idx_row,current+1,n,row_hint,tmp,row_pattern);
            pattern.push_back('0');
        }
    }

}



void generate_col(int idx_col,int current,int n,vector<int>&col_hint, string pattern, vector< vector< string > > &col_pattern)//according to the row hint to generate_row
{   
    if (col_hint.empty())
    {
        // cout<<"empty"<<endl;
        string tmp="";
        for(int i =0;i<n;i++)
        {
            tmp.push_back('0');
        }
        col_pattern[idx_col].push_back(tmp);
        for(int i=0;i<n;i++)
        {
            guess_board[i][idx_col]='0';
        }
        return;
    }
    else if(current>=col_hint.size())
    {
        // cout<<"current<size"<<endl;
        int space=n-pattern.size();
        if(space<0)
        {
            pattern.pop_back();
        }
        else
        {
            for(int i=0;i<space;i++)
            {
                pattern.push_back('0');
            }
        }
        col_pattern[idx_col].push_back(pattern);
        return;
    }
    else
    {   
        // cout<<"current:"<<current<<endl;
        int sum=0;//sum: # need to be 1
        for(int i=current;i<col_hint.size();i++)
        {
            sum+=col_hint[i];
        }
        int space=n-pattern.size()-sum+col_hint.size()-current-1;//row_hint-current-1 : necessary space between left pattern 
        // cout<<"space:"<<space<<endl;
        if(space<0){return;}
        for(int i=0;i<=space;i++)
        {
            string tmp=pattern;
            for(int j=0;j<col_hint[current];j++)
            {
                tmp.push_back('1');
            }
            tmp.push_back('0');
            generate_col(idx_col,current+1,n,col_hint,tmp,col_pattern);
            pattern.push_back('0');
        }
    }

}



bool guess(int &n,vector< vector<string > > &row_pattern, vector< vector<string > > &col_pattern)
{
    bool change=false;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            char tmp=row_pattern[i][0][j];
            bool diff=false;
            for(int k=1;k<row_pattern[i].size();k++)
            {
                if(tmp!=row_pattern[i][k][j])
                {
                    diff=true;
                    break;
                }
            }
            if(!diff)
            {
                guess_board[i][j]=tmp;
                diff=false;
            }

            tmp=col_pattern[i][0][j];
            for(int k=1;k<col_pattern[i].size();k++)
            {
                if(tmp!=col_pattern[i][k][j])
                {
                    diff=true;
                    break;
                }
            }
            if(!diff)
            {
                guess_board[j][i]=tmp;
                diff=false;
            }
        }
    }

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(guess_board[i][j]!='*')
            {
                if(eliminate(i,j,guess_board[i][j],row_pattern,col_pattern))
                {
                    change=true;
                }
                
            }
        }
    }
    return change;
}




int main(int argc, char * argv[])
{
    clock_t t1,t2;
    double totaltime=0.0;;
   
    t1=clock();
    ifstream file;
    fstream output;
    
    string line;
    if(argc==1)
    {
        file.open("15_all.txt",ios::in);
        output.open("result.txt",ios::out);
    }
    else
    {
        file.open(argv[1],ios::in);
        output.open(argv[2],ios::out);
    }

    
    int line_num=30;

    
    
    file.clear();
    file.seekg(0,ios::beg);
    int n=line_num/2;
    
    string init;
    for(int i=0;i<n;i++)
    {
        init.push_back('*');
    }
    for(int i=0;i<n;i++)
    {
        guess_board.push_back(init);
    }
    
    vector< vector< vector<int> > > hint_cols;
    vector< vector< vector<int> > > hint_rows;

    vector<vector<string > > row_pattern;
    vector<vector<string > > col_pattern;
    int count=0;
    file.seekg(0,ios::beg);
    while(getline(file,line))
    {
        if(line[0]=='$')
        {
            
            
            vector< vector<int> > col;
            vector< vector<int> > row; 
            for(int i=0;i<n;i++)
            {
                vector<int>tmp;
                getline(file,line);
                stringstream ss(line);
                int num;
                while(ss>>num)
                {
                    
                    tmp.push_back(num);
                }
                col.push_back(tmp);
            }

            for(int i=0;i<n;i++)
            {
                vector<int>tmp;
                getline(file,line);
                stringstream ss(line);
                int num;
                while(ss>>num)
                {
                    tmp.push_back(num);
                   
                }
                row.push_back(tmp);
            }
            count++;
            vector<string> tmp;
            string tmpv="";
            // cout<<'w'<<endl;
            for (int i=0;i<n;i++)
            {
                
                row_pattern.push_back(tmp);
                col_pattern.push_back(tmp);
                generate_row(i,0,n,row[i],tmpv,row_pattern);
                generate_col(i,0,n,col[i],tmpv,col_pattern);
                // cout<<"row_patternsize:"<<row_pattern[i].size()<<endl;
                // cout<<"col_patternsize:"<<col_pattern[i].size()<<endl;
            }

            
            for(int i=0;i<n;i++)
            {
                if(!guess(n,row_pattern,col_pattern))
                {
                    break;
                }
            }
            
            long long int complex_row=1;
            long long int complex_col=1;
            for(int i=0;i<n;i++)
            {
                
                complex_row*=row_pattern[i].size();
                complex_col*=col_pattern[i].size();
                if(complex_row<0||complex_col<0)
                {
                    break;
                }
            }

            
            cout<<"prob:"<<count<<endl;
            vector<string> board;
            
            
            string result="";
            DFS(0,n,col,board,row_pattern,col_pattern,result);
            output<<"$"<<count<<endl;
            output<<result;


            
            // cout<<"board"<<endl;
            //clear pattern for next problem
            for(int i=0;i<n;i++)
            {
                row_pattern[i].clear();
                col_pattern[i].clear();
                guess_board[i]=init;
            }
            
            t2=clock();
            cout<<"subtime:"<<(t2-t1)/(double)(CLOCKS_PER_SEC)<<endl;
            totaltime+=(t2-t1)/(double)(CLOCKS_PER_SEC);
            t1=t2;
            // break;
        }    
        // cout<<ss;
    }
    cout<<endl<<"time:"<<totaltime<<endl;
    file.close();  
    output.close();
}