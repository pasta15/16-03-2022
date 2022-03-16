import os
from flask import Flask, render_template,redirect, url_for,request
import re
global cation
global anion
global enthalpy_type
global message
import sqlite3
import csv


app = Flask(__name__)




def tuple2list(aList):
  newlist=[]
  for i in aList:
    newlist.append(i[0])
  return newlist

def setup():
    global valid_cations
    global trans_metal
    global valid_elements
    global valid_anions
    con = sqlite3.connect("data.db")
    cur = con.cursor()
    query_result = cur.execute("SELECT * FROM cation")
    valid_cations=dict(query_result)

    query_result = cur.execute("SELECT * FROM anion")
    valid_anions=dict(query_result)

    query_result = cur.execute("SELECT * FROM TransMetals")
    temp_transmetal=list(query_result)
    trans_metal=tuple2list(temp_transmetal)

    query_result = cur.execute("SELECT * FROM element")
    valid_elements=dict(query_result)




    con.commit()

def validation(usalt):
    global cation
    global anion
    cation=''
    anion=''
    usalt=usalt.replace(' ' ,'') #remove spaces
    if ord(usalt[0]) not in range(65,91):
        return False
    else:
        usalt1=usalt[1:len(usalt)] #this will remove the first character as we know it will be a capital letter
        position = re.search("[A-Z]",usalt1) #this will find the position of the next capital letter
        if position==None:
            return False
        else:
            cation=usalt[0:position.start()+1]
            anion=usalt[position.start()+1:len(usalt)]
            #cation and anion will store the values of the first and last half of the salt
            return True

def validcation():
    global ccharge
    global cation1
    cation1=cation
    #valid_cations={'Al':3,'H':1,'Li':1,'Na':1,'K':1,'Be':2,'Mg':2,'Ca':2,'Fe':0,'Cu':0,'Zn':0}
    #trans_metal=['Fe','Cu','Zn']
    global valid_cations
    global trans_metal
    global IsTransmetal
    IsTransmetal=False
    global cbalance
    cbalance=1
    cation1 = re.sub('\d', '', cation1[::-1],1)  #this will create a copy of the metal the last number (if there is one)
    cbalance=re.findall('\d',cation) #this will store the last number in the salt as they will be needed to check the salt is balanced. E.g in Na2O it will store the 2 for later. cbalance stands for cation balancing number.
    if cbalance==[]: #i.e. if it was empty/the cation has no numbers
        cbalance=1
    else:
        cbalance=cbalance[0]
    #Any extra numbers will cause the validation to fail.
    cation1 = re.sub('\(', '', cation1)  #if there is a polyatomic nonmetal then there might be an opening bracket so this line will remove it. e.g for Mg(OH)2 cation1 would be Mg(
    cation1=cation1[::-1]
    if cation1 in valid_cations:
        if cation1 in trans_metal:
            IsTransmetal=True
            return True
        else:
           ccharge=valid_cations[cation1]
           return True
    
    else:
        return False


global valid_elements
#valid_elements={'N':3,'O':2,'S':2,'F':1,'Cl':1,'Br':1,'I':1}
global valid_anions
#valid_anions={'N':3,'O':2,'S':2,'F':1,'Cl':1,'Br':1,'I':1,'SO4':2,'OH':1,'NO3':1,'PO4':3}

######
def remove_brackets():
    global cbalance
    global anion
    global anion1
    global ccharge
    last=len(anion)
    if anion[last-1]!= ')':#this checks whether the last character is a bracket because if it isn�t then the user has not provided enough information to construct a Born-Haber cycle.
        return False
    position = re.search('\(',anion)
    count=0
    for i in range(position.start()+1,last-1) :
        if anion[i]== 'I':
            count=count+1
        else:
           return False
    #this loop will count the number of Is. So for (III) it will work out that the charge will be 3

    ccharge=count
    anion1=re.sub('I','',anion) #this will remove the brackets and roman numerals so that the non-metal can be validated
    anion1=re.sub('\(','',anion1[::-1],1) #For these two lines the program will remove the first set of () starting at the end. I have done this incase there are brackets that have been used to balance a polyatomic anion.
    anion1=re.sub('\)','',anion1[::-1],1)

######

#####
def validanion():
    global enthalpy_type
    global anion
    global acharge
    global abalance
    global anion1
    abalance=1
    if IsTransmetal==True:
        if remove_brackets()==False:
            return False
    else:
        anion1=anion
    if ')' in anion1:
        if anion1[len(anion1)-1] !=')':
            anion1=list(anion1)
            abalance=anion1.pop()
            anion1=''.join(anion1)
    else:
        if (ord(anion1[len(anion1)-1]) in range(49,58)):  #ie if the anion is an element and not a poly atomic one. and it does have an extra number
            if anion1 not in valid_anions:  
                anion1=list(anion1)
                abalance=anion1.pop()
                anion1=''.join(anion1)

     
    anion1  = re.sub("\)", "", anion1 ,1) #removes the bracket  
    if enthalpy_type=='solution': #this if statement is needed as there is a different list of valid non-metals for enthalpy change of solution and formation Born-Haber cycles.
        if anion1 in valid_anions:
           acharge=valid_anions[anion1]
           return True
        else:            
            return False
    else:
       anion1= re.sub('\d','',anion1)  #this will remove any numbers left
       abalance=re.findall('\d',anion)
    if abalance==[]: #i.e. if it was empty/the cation has no numbers
       abalance=1
    else:
        abalance=abalance[0]
    #this will save the balancing number (if any)
    if anion1 in valid_elements:
        acharge=valid_elements[anion1]
        return True
    else:
       return False
####

####
def balanced():
    global abalance
    global cbalance
    abalance=int(abalance)
    cbalance=int(cbalance)
    negative=(acharge*abalance)

    positive=(ccharge*cbalance)
    if (negative)==(positive):
        return True
    else:
        return False

###
def enthalpychange():
    global salt_type
    global valid_enthalpies
    global valid_enthalpies_names
    global enthalpy_change
    global enthalpy_type
    global location
    global noof
    

    con = sqlite3.connect("data.db")
    cur = con.cursor()
    if enthalpy_type=='solution':
        salt_type='sol'
        query_result_id = cur.execute("SELECT change FROM changes WHERE type='%s' " %'sol')
        valid_enthalpies=tuple2list(query_result_id)
        query_result_full = cur.execute("SELECT fullname FROM changes WHERE type='%s' " %'sol') #new line
        valid_enthalpies_names=tuple2list(query_result_full)
    else:
        salt_type=str(cbalance)+str(ccharge)+str(abalance)+str(acharge)
        query_result_id = cur.execute("SELECT change FROM changes WHERE type='%s' OR type='0'" %salt_type)
        valid_enthalpies=tuple2list(query_result_id)
        query_result_full = cur.execute("SELECT fullname FROM changes WHERE type='%s' OR type='0'" %salt_type) #new line
        valid_enthalpies_names=tuple2list(query_result_full)

    #valid_enthalpies=tuple2list(list(query_result_id))
    query_result_noof=cur.execute("SELECT noofsteps FROM salts WHERE type='%s'" %salt_type)
    noof=list(query_result_noof)
    noof=tuple2list(noof)
    noof=noof[0]
    

    for i in range(0,len(valid_enthalpies_names)):
        valid_enthalpies_names[i]=valid_enthalpies_names[i].replace('cation',cation1)
        valid_enthalpies_names[i]=valid_enthalpies_names[i].replace('anion',anion1)
        valid_enthalpies_names[i]=valid_enthalpies_names[i].replace('salt',salt1)


    con.close()
    
    if enthalpy_change in valid_enthalpies:
        location=valid_enthalpies.index(enthalpy_change)
        return True
    else:
        return False
    


###
####
def run_validation(asalt):
    global message
    message='an invalid salt.\\nThe most likely reason as to why is because '
    if validation(asalt)==False:
        message=message + 'the format not understandable'
        return False
    if validcation()==False:
        message=message + 'the cation is unrecognisable'
        return False
    if validanion()==False:
        message=message + 'the anion is unrecognisable' 
        return False
    if balanced()==False:
        message=message + 'the salt is unbalanced'
        return False
    if enthalpychange()==False:
        message=message + 'the salt or the selected type of enthalpy cycle does not have that enthalpy change'
        return False
   
    else:
        return True




@app.route("/",methods=["POST","GET"])  # this sets the route to this page
def home():
    setup()
    global message
    global enthalpy_type
    global salt1
    global enthalpy_change
    if request.method == "POST":
        salt1 = request.form["salt"]
        enthalpy_type= request.form["type"]
        enthalpy_change=request.form["change"]
        if run_validation(salt1)==False:
            return render_template("intro.html",content=message,need='true')
        else:

            return redirect(url_for('mid'))
            ##return render_template("mid.html",len=len(valid_enthalpies),changes=valid_enthalpies)
        #return render_template("intro.html",content=message,need='true')
        #if run_validation(salt1)==False:
        #    test=run_validation(salt1)
       #     test1=message
       #     return render_template("intro.html",content=message,need='true')
      #      
       # else:
        #    return render_template("intro.html",content=message,need='true')
        #return render_template("intro.html",content=" alert('Hello! I am an alert box!!');")
    else:
        return render_template("intro.html",need='')
        #
        #return render_template("intro.html",content=" alert('Hello! I am an alert box!!');")

	#return render_template("intro.html") 



def sol_math():
    global total_solution
    global equation
    global values
    values=[]
    divider=1
    a_mult=''
    c_mult=''
    if abalance!=1:
        a_mult=str(abalance) + 'x'

    if cbalance!=1:
       c_mult=str(cbalance) + 'x'
    values.append(equation[0])
    values.append(c_mult + str(equation[1]))
    values.append(a_mult + str(equation[2]))
    values.append(equation[3])
    if equation[1]=='x':
        divider=cbalance
    else:
        equation[1]=equation[1]*cbalance

    if equation[2]=='x':
        divider=abalance
    else:
        equation[2]=equation[2]*abalance
    #####
    sum1=0
    answer1=0
    if equation.index('X')==0:
        sum1=equation[1]+equation[2]
        answer1=sum1-equation[3]

    if equation.index('X')==3:
        sum1=equation[1]+equation[2]
        answer1=sum1-equation[0]

    if equation.index('X')==1:
        sum1=equation[0]+equation[3]
        answer1=sum1-equation[2]
        answer1=answer1/divider

    if equation.index('X')==2:
        sum1=equation[0]+equation[3]
        answer1=sum1-equation[1]
        answer1=answer1/divider

    #print(answer1)

    #he first line in each case statement adds up the numbers on the other side of the ‘x’. The second line solves for ‘x’. 
    #I used a case statement as it is the easiest way to have different lines that run depending on the value of one variable. In addition it is more understandable than several if else statements.

    #Procedure to output the calculations for solution born haber cycles
   # setup=[]
    #print('s ',setup)
    setup=(str(equation[0]) + ' + ' + str(equation[3]) + ' = ' + c_mult + str(equation[1]) + ' + ' + a_mult + str(equation[2]))
    #print(setup)
    total_solution=setup + '\n'
    if equation.index('X')==0:
     total_solution=total_solution + ' X ' + str(equation[3]) + ' =' + str(sum1) + '\n'
    if equation.index('X')==1:
        total_solution=total_solution + 'X' + str(equation[2]) + '=' + str(sum1) +'\n'
    if equation.index('X')==2:
        total_solution=total_solution + 'X' + str(equation[1]) + '=' + str(sum1) +'\n'
    if equation.index('X')==3:
        total_solution=total_solution + 'X' + str(equation[0]) + '=' + str(sum1) + '\n'

    if divider!=1:
      total_solution=total_solution + str(answer1*divider) + ' = ' + 'X' + '\n'
      total_solution=total_solution.replace(str(divider)+ ' x X','X',1)
      total_solution=total_solution.replace('X',str(divider)+ 'X')

    total_solution=total_solution + 'x' + '=' +str(answer1)
    total_solution=total_solution.replace('\n','<br>')
          
#print('x' + '=' +str(answer1))




def form_math():
    global equation
    global total_solution
    global values
    global cbalance
    divider=1
    values=[]
    setup=''
    values.append(equation[0])
    a_mult=''
    c_mult=''
    temp=equation
    setup=str(equation[0]) + ' = '
    if abalance!=1:
        a_mult=str(abalance) + 'x'

    if cbalance!=1:
       c_mult=str(cbalance) + 'x'

    values.append(c_mult+str(equation[1]))
    values.append(a_mult+str(equation[2]))
        #the first value (starting from the bottom left) never needs to be multiplied.
    if equation[1]=='X':
    #the second value needs to be multiplied by the nonmetal’s balancing number
      divider=cbalance
      setup=setup + c_mult+ str(equation[1]) + ' + ' 
    else:
      setup=setup + c_mult+ str(equation[1]) + ' + ' 
      equation[1]=equation[1] * cbalance

    if equation[2]=='X':
    #the second value needs to be multiplied by the nonmetal’s balancing number
      divider=abalance
      setup=setup + a_mult+ str(equation[2]) + ' + ' 
    else:
      setup=setup + a_mult+ str(equation[2]) + ' + ' 
      equation[2]=equation[2] * abalance

    temp1=equation
    for i in range(3,3+ccharge): ##The loop will start at the 4th enthalpy change (which is always first ionization energy) and do extra steps based on the charge of the cation because the charge dictates how many ionization changes there are
      if equation[i]=='X':
        divider=cbalance
      else:
        setup=setup + c_mult+ str(equation[i]) + ' + ' 
        equation[i]=equation[i] * cbalance
      values.append(c_mult+str(equation[i]))

    for i in range (3+ccharge,3+ccharge+acharge):
      #The loop will start at the enthalpy change after the ionization energies and do extra steps based on the charge of the anion because the charge dictates how many electron affinity changes there are
      if equation[i]=='X':  
        divider=abalance
      else:
         setup=setup + a_mult+ str(equation[i]) + ' + ' 
         equation[i]=equation[i] * abalance
      values.append(a_mult+str(equation[i]))
    temp2=equation
    sum=0
    answer=''
    for i in range(1,len(equation)): #its from 1 as we want to skip the first one
      if equation[i]!='X':
        sum=sum+equation[i] #this loop will add up all the numbers on the second route

    if equation[0]=='X':  		 #this if statement checks if we have to rearrange the equation
      answer=sum	 
    else :
        answer=(equation[0]-sum)/divider

    answer=int(answer)


    total_solution=''
    setup=setup+str(equation[len(equation)-1])
    setup=str(equation[0]) + '=' 
    for i in range(1,len(equation)-1): #the  -1 is  because the list starts at zero
        setup=setup  + str(equation[i]) + '+'
  

    #temp=len(equation)-1
    setup=setup+str(equation[len(equation)-1])
   # print(setup)
    total_solution=total_solution+setup +'\n'

    if equation[0] !='X':
        total_solution=total_solution+ str(equation[0]) + '='  + 'X' + ' + '+ str(sum) + '\n'
        total_solution= total_solution+str(equation[0]) + ' - '  + str(sum) + ' = ' + 'X' +'\n'
        if divider==1:
            total_solution=total_solution+'X' + '=' + str(answer)
        else:
            total_solution=total_solution + str(answer*divider) + ' = ' + 'X' + '\n'
            total_solution=total_solution.replace(str(divider)+ ' x X','X',1)
          #the above line will remove the section the X being multiplied so it can be replaced with divider+X
      
        total_solution=total_solution.replace('X',str(divider)+ 'X')
        total_solution=total_solution+'X' + '=' + str(answer)
    else:
        total_solution=total_solution+'X' + '=' + str(answer)
        

    total_solution=total_solution.replace('\n', '<br>')

###
def math():
    global noof
    global location
    global equation
    noofsteps= noof#total#the number of steps in the Born-Haber cycle
 
    equation=[] #blank list
 
    count=1
    poslist=0
 
    while count<=noofsteps: #this loop will put the values from the user and an ‘x’ in the right position so that the later procedures can solve for x.
      if count==location+1:
        equation.append('X')
      else:
        equation.append(data[poslist])
        poslist=poslist+1 
        #print('poslitst',poslist)#this is used to make sure no item in data is skipped
      count=count+1
    equation1=equation.copy()
    
     # print(count)
  
#    print(equation)
    if enthalpy_type=='solution':
        sol_math()
    else:
        form_math()





####
   

###



def valid_data():
    global data
    global message1

    for i in range(0,len(data)):
        try:
            data[i]=int(data[i])
            if data[i] not in range (-7500,7500):
                message1='are out of range'
                return False
        except:
            message1='contain non-numerical digits'
            return False
    return True



@app.route("/data",methods=["POST","GET"])  # this sets the route to this page
def mid():
    global data
    global message1
    global total_solution
    global cbalance
    global values
    global diatomics
    diatomics=['N','O','Br','F','Cl','I']
    if request.method == "POST":
        data=[]
        for i in request.form:
            data.append(request.form[i])
        if valid_data()==True:
            math()
            template=salt_type+".html"
            atempbalance=''
            cbalance=''
            if abalance!=1:
                atempbalance=abalance
            if cbalance!=1:
                ctempbalance=cbalance
            if anion1 in diatomics:
                tempanion=anion1+'2'
                temp_atempbalance=abalance/2
            else:
                tempanion=anion1
                temp_atempbalance=abalance


            return render_template(template,cation=cation1,anion=anion1,salt=salt1,abalance=atempbalance,cbalance=ctempbalance,values=values,answer=total_solution,tempanion=tempanion,tempbalance=temp_atempbalance)
        else:
           return render_template("mid.html",len=len(valid_enthalpies),changes=valid_enthalpies_names,pos=location,content1=message1,need1='true',cation='na')

      #  if valid_data==False:
       #     return #error
       # else:
        #    return render_template("base.html")


    else:
        return render_template("mid.html",len=len(valid_enthalpies),changes=valid_enthalpies_names,pos=location,need='')


if __name__ == "__main__":
    #app.run(debug=True)
    app.run()


 #def x():
 #return render_template("mid.html",len=len(valid_enthalpies),changes=valid_enthalpies_names,pos=location),content1=message1,need1='true')
