## author haowang 
import numpy as np
check_pair_dict={1:31,2:32,3:34,4:33,5:36,6:35,7:25,8:26,9:28,10:27,11:30,12:29,13:37,14:38,15:40,16:39,17:45,18:46,19:48,20:47,21:41,22:42,23:44,24:43,25:7,26:8,27:10,28:9,29:12,30:11,31:1,32:2,33:4,34:3,35:6,36:5,37:13,38:14,39:16,40:15,41:21,42:22,43:24,44:23,45:17,46:18,47:20,48:19}

def nearest_hopping():
    with open(r"C:\Users\wangh\Desktop\WF1_hr.dat","r") as f:
        lines=f.readlines()
        wan_num=int(lines[1]); wsc_num=int(lines[2])
        ski_row_num=int(np.ceil(wsc_num/15.0))
        with open(r"C:\Users\wangh\Desktop\WF1_hr_1_032_1st2nd.dat","a") as f1:
            for i in range(ski_row_num+3):
                f1.write(lines[i])
        for i in range(wan_num**2*wsc_num):
            if abs(float(lines[3+ski_row_num+i].split()[0])) > 2  or abs(float(lines[3+ski_row_num+i].split()[1])) > 2 :
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[5]),"0.000000")
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[6]),"0.000000")

#### remove_Sr_p
            if abs(float(lines[3+ski_row_num+i].split()[3])) >= 14 and abs(float(lines[3+ski_row_num+i].split()[3])) <= 16 :
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[5]),"0.000000")
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[6]),"0.000000")

            if abs(float(lines[3+ski_row_num+i].split()[4])) >= 14 and abs(float(lines[3+ski_row_num+i].split()[4])) <= 16 :
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[5]),"0.000000")
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[6]),"0.000000")

            if abs(float(lines[3+ski_row_num+i].split()[3])) >= 38 and abs(float(lines[3+ski_row_num+i].split()[3])) <= 40 :
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[5]),"0.000000")
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[6]),"0.000000")

            if abs(float(lines[3+ski_row_num+i].split()[4])) >= 38 and abs(float(lines[3+ski_row_num+i].split()[4])) <= 40 :
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[5]),"0.000000")
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[6]),"0.000000")
#### remove_Sr_p

#####0.01
            if abs(float(lines[3+ski_row_num+i].split()[5])) <= 0.032 :
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[5]),"0.000000")

            if abs(float(lines[3+ski_row_num+i].split()[6])) <= 0.032 :
                lines[3+ski_row_num+i]=lines[3+ski_row_num+i].replace(str(lines[3+ski_row_num+i].split()[6]),"0.000000")

            with open(r"C:\Users\wangh\Desktop\WF1_hr_1_032_1st2nd.dat","a") as f1:
                f1.write(str(lines[3+ski_row_num+i]))

#nearest_hopping()


def inversion_operation():
    with open(r"C:\Users\wangh\Desktop\WF1_hr_1_032_1st2nd.dat","r") as f:
        lines=f.readlines()
        wan_num=int(lines[1]); wsc_num=int(lines[2])
        ski_row_num=int(np.ceil(wsc_num/15.0))
        Hmn_R = []
        for i in range(wan_num**2*wsc_num): 
            for j in range(7):
                Hmn_R.append(float(lines[3+ski_row_num+i].split()[j]))

        Hmn_R_np = np.array(Hmn_R)
        Hmn_R_np = Hmn_R_np.reshape((wan_num**2*wsc_num,7))

        for i in range(wsc_num):
            for j in range(wsc_num):
                if Hmn_R_np[wan_num**2*i][0] == -1 * Hmn_R_np[wan_num**2*j][0] and Hmn_R_np[wan_num**2*i][1] == -1 * Hmn_R_np[wan_num**2*j][1] and abs(Hmn_R_np[wan_num**2*i][0]) <= 2 and abs(Hmn_R_np[wan_num**2*j][0]) <= 2 and abs(Hmn_R_np[wan_num**2*i][1]) <= 2 and abs(Hmn_R_np[wan_num**2*j][1]) <= 2:
                    # print(Hmn_R_np[wan_num**2*i][0],Hmn_R_np[wan_num**2*i][1])
                    # print(Hmn_R_np[wan_num**2*j][0],Hmn_R_np[wan_num**2*j][1])
                    for k in range(wan_num**2):
                        for l in range(wan_num**2):
                            if Hmn_R_np[wan_num**2*i+k][3] == Hmn_R_np[wan_num**2*j+l][4] and Hmn_R_np[wan_num**2*i+k][4] ==Hmn_R_np[wan_num**2*j+l][3]:
                                average_real = (Hmn_R_np[wan_num**2*i+k][5]+Hmn_R_np[wan_num**2*j+l][5])/2
                                Hmn_R_np[wan_num**2*i+k][5] = average_real
                                Hmn_R_np[wan_num**2*j+l][5] = average_real
                                average_image = (abs(float(Hmn_R_np[wan_num**2*i+k][6])) + abs(float(Hmn_R_np[wan_num**2*j+l][6])))/2
                                if Hmn_R_np[wan_num**2*i+k][6] > 0:
                                    Hmn_R_np[wan_num**2*i+k][6] = average_image
                                else:
                                    Hmn_R_np[wan_num**2*i+k][6]= -1 * average_image
                                if Hmn_R_np[wan_num**2*j+l][6] > 0:
                                    Hmn_R_np[wan_num**2*j+l][6] = average_image
                                else:
                                    Hmn_R_np[wan_num**2*j+l][6] = -1 * average_image

        np.savetxt(r"C:\Users\wangh\Desktop\WF1_hr_1_032_1st2nd_inv.dat",Hmn_R_np,fmt=['%5d']*5+['%15f']*2)
inversion_operation()

def inversion2_operation():
    with open(r"C:\Users\wangh\Desktop\WF1_hr_2.dat","r") as f:
        lines=f.readlines()
        wan_num=int(lines[1]); wsc_num=int(lines[2])
        ski_row_num=int(np.ceil(wsc_num/15.0))
        Hmn_R = []
        for i in range(wan_num**2*wsc_num): 
            for j in range(7):
                Hmn_R.append(float(lines[3+ski_row_num+i].split()[j]))

        Hmn_R_np = np.array(Hmn_R)
        Hmn_R_np = Hmn_R_np.reshape((wan_num**2*wsc_num,7))

        for i in range(wsc_num):
            for j in range(wsc_num):
                if Hmn_R_np[wan_num**2*i][0] == 0 and Hmn_R_np[wan_num**2*j][0] == 0 and Hmn_R_np[wan_num**2*i][1] == 0 and Hmn_R_np[wan_num**2*j][1] == 0:
                    for k in range(wan_num**2):
                        if Hmn_R_np[wan_num**2*i + k][3] == Hmn_R_np[wan_num**2*i + k][4]:
                            real_1 = Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][5]
                            real_2 = Hmn_R_np[wan_num**2*i + k][5]
                            real_aver = (real_1+real_2)/2
                            Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]] - 1 )* 48  + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][5] = real_aver
                            Hmn_R_np[wan_num**2*i + k][5] =real_aver
                            image_1 = Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][6]
                            image_2 = Hmn_R_np[wan_num**2*i + k][6]
                            image_aver = (image_1 + image_2)/2
                            Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] - 1][6] = image_aver
                            Hmn_R_np[wan_num**2*i + k][6] = image_aver
        np.savetxt(r"C:\Users\wangh\Desktop\WF1_hr_3.dat",Hmn_R_np,fmt=['%5d']*5+['%15f']*2)
#inversion2_operation()

#check_pair_dict={1:31,2:2,3:3,4:4,5:5,6:6,7:25,8:8,9:9,10:10,11:11,12:12,13:37,14:38,15:40,16:39,17:45,18:46,19:48,20:47,21:41,22:42,23:44,24:43,25:7,26:26,27:27,28:28,29:29,30:30,31:1,32:32,33:33,34:34,35:35,36:36,37:13,38:14,39:16,40:15,41:21,42:22,43:24,44:23,45:17,46:18,47:20,48:19}
def mirror_operation():
    with open(r"C:\Users\wangh\Desktop\WF1_hr_2.dat","r") as f:
        lines=f.readlines()
        wan_num=int(lines[1]); wsc_num=int(lines[2])
        ski_row_num=int(np.ceil(wsc_num/15.0))
        Hmn_R = []
        for i in range(wan_num**2*wsc_num): 
            for j in range(7):
                Hmn_R.append(float(lines[3+ski_row_num+i].split()[j]))

        Hmn_R_np = np.array(Hmn_R)
        Hmn_R_np = Hmn_R_np.reshape((wan_num**2*wsc_num,7))

    for i in range(wsc_num):
        for j in range(wsc_num):
            if Hmn_R_np[wan_num**2*i][0] == -1 * Hmn_R_np[wan_num**2*j][1] and Hmn_R_np[wan_num**2*i][1] == -1 * Hmn_R_np[wan_num**2*j][0] and abs(Hmn_R_np[wan_num**2*i][0]) <= 2 and abs(Hmn_R_np[wan_num**2*j][0]) <= 2 and abs(Hmn_R_np[wan_num**2*i][1]) <= 2 and abs(Hmn_R_np[wan_num**2*j][1]) <= 2:
                # print(Hmn_R_np[wan_num**2*i][0],Hmn_R_np[wan_num**2*i][1])
                # print(Hmn_R_np[wan_num**2*j][0],Hmn_R_np[wan_num**2*j][1])
                for k in range(wan_num**2):                  
                    real_1 = Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][5]
                    real_2 = Hmn_R_np[wan_num**2*i + k][5]
                    real_aver = (real_1+real_2)/2
                    Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]] - 1 )* 48  + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][5] = real_aver
                    Hmn_R_np[wan_num**2*i + k][5] =real_aver
                    image_1 = Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][6]
                    image_2 = Hmn_R_np[wan_num**2*i + k][6]
                    image_aver = (image_1 + image_2)/2
                    Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] - 1][6] = image_aver
                    Hmn_R_np[wan_num**2*i + k][6] = image_aver

        np.savetxt(r"C:\Users\wangh\Desktop\results.txt",Hmn_R_np,fmt=['%5d']*5+['%15f']*2)

#mirror_operation()

def c2_operation():
    with open(r"C:\Users\wangh\Desktop\WF1_hr_1.dat","r") as f:
        lines=f.readlines()
        wan_num=int(lines[1]); wsc_num=int(lines[2])
        ski_row_num=int(np.ceil(wsc_num/15.0))
        Hmn_R = []
        for i in range(wan_num**2*wsc_num): 
            for j in range(7):
                Hmn_R.append(float(lines[3+ski_row_num+i].split()[j]))

        Hmn_R_np = np.array(Hmn_R)
        Hmn_R_np = Hmn_R_np.reshape((wan_num**2*wsc_num,7))

    for i in range(wsc_num):
        for j in range(wsc_num):
            if Hmn_R_np[wan_num**2*i][0] == -1 * Hmn_R_np[wan_num**2*j][1] and Hmn_R_np[wan_num**2*i][1] == -1 * Hmn_R_np[wan_num**2*j][0] and abs(Hmn_R_np[wan_num**2*i][0]) <= 2 and abs(Hmn_R_np[wan_num**2*j][0]) <= 2 and abs(Hmn_R_np[wan_num**2*i][1]) <= 2 and abs(Hmn_R_np[wan_num**2*j][1]) <= 2:
                # print(Hmn_R_np[wan_num**2*i][0],Hmn_R_np[wan_num**2*i][1])
                # print(Hmn_R_np[wan_num**2*j][0],Hmn_R_np[wan_num**2*j][1])
                for k in range(wan_num**2):
                    # print(Hmn_R_np[wan_num**2*j + check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]] * 48 -1 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]]][0])          
                    # print(Hmn_R_np[wan_num**2*j + check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]] * 48 -1 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]]][1])
                    # print(Hmn_R_np[wan_num**2*j + check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]] * 48 -1 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]]][3])
                    # print(Hmn_R_np[wan_num**2*i + k][4])                  
                    
                    real_1 = Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][5]
                    real_2 = Hmn_R_np[wan_num**2*i + k][5]
                    real_aver = (real_1+real_2)/2
                    Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]] - 1 )* 48  + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][5] = real_aver
                    Hmn_R_np[wan_num**2*i + k][5] =real_aver
                    image_1 = Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] -1 ][6]
                    image_2 = Hmn_R_np[wan_num**2*i + k][6]
                    image_aver = (image_1 + image_2)/2
                    Hmn_R_np[wan_num**2*j + (check_pair_dict[Hmn_R_np[wan_num**2*i+k][3]]-1) * 48 + check_pair_dict[Hmn_R_np[wan_num**2*i+k][4]] - 1][6] = image_aver
                    Hmn_R_np[wan_num**2*i + k][6] = image_aver

        np.savetxt(r"C:\Users\wangh\Desktop\results.txt",Hmn_R_np,fmt=['%5d']*5+['%15f']*2)

#c2_operation()