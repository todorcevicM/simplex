from xml.etree.ElementPath import find
import numpy as np


def simplex(A, b, c, max = True):
    m, n = A.shape
    table = np.zeros((m + 1, m + 3))

    # prvi korak uvek isti, uvek kao pocetne bazne promenljive biramo dodatne
    # koliko treba dodati promenljivih odredjuje koliko ce biti baznih promenljivih
    # broj ogranicenja je m
    # broj promenljivih je n
    # ali zapravo broj promenljivih raste na n + m
    # promenljivih u c mora biti isto koliko i u A
    if (n != c.shape[0]):
        print("Broj promenljivih niji isti!")
        return
    if (m != b.shape[0]):
        print("Nije dobro zadato ogranicenje!")
        return
    
    # govori koje su tacno bazne promenljive
    # 0 nije bazna, 1 jeste bazna
    base = np.zeros(n, dtype=int)
    # dodajem promenljive koje nisu u originalnim ogranicenjima i samim tim pretvaram nejednacine u jednacine
    # te promenljive su uvek na kraju i 1 oznacava da su bazne
    for _ in range(m):
        base = np.append(base, 1)
        # takodje potrebno je dodati te novododate promenljivi u originalni c-vekor jer i one sada uticu na max problem
        # posto nisu bile u origalnom max problemu, dodaju se 0 jer te promenljive uticu 0
        c = np.append(c, 0)

    B = np.identity(m)    
    A = np.append(A, B, axis=1)

    # TODO: ovde izmeniti ako se radi primer gde one bazne promenljive nisu na dodatne promenljive
    # base = np.array([1, 1, 0, 1, 0])
    # B = np.matrix([[-1 / 2, -1, 0], [-2, -2, 1], [-1, -4, 0]])

    # na pocetku mogu kreirati B i B_inv matricu
    # sada moze da se krene u obradu
    # posto cu update-ovati base-vektor na kraju iteracije
    B_inv = np.linalg.inv(B)
    while True:
        # tabela
        # nulta kolona su bazne promenljive
        j = 0
        k = 0
        for i in base:
            j += 1
            if i == 1:
                k += 1
                table[k: k + 1, 0] = j
        
        print(table)
        print(B_inv)
        # sledecih n kolona je B_inv
        table[1:, 1:-2] = B_inv
        # sledeca slobodna kolona je slobodna kolona
        # slobodna kolona
        b_B = B_inv @ b
        table[1:, -2] = b_B.transpose()

        # nulti red, od prve kolone do dve kolone od nazad je omega
        # za omega mi je potrebno c_B koji sadrzi koeficijente baznih promenljivih

        c_B = [i for (index, i) in enumerate(c) if base[index] == 1]
        print(f"c_B vektor: {c_B}")
        w = c_B @ B_inv
        print(f"w vektor: {w}")
        table[0, 1:-2] = w

        # nakon sto imamo w, mozemo izracunati i vrednost u nultom redu i slobodnoj koloni
        # ta vrednost se dobija kao: w @ b_B
        table[0, -2] = c_B @ b_B

        # sledece se racunaja pivot kolona
        # neophodno je naci najmanji element u nizu koji se dobija ubacivanjem primenom operacije na koeficijente pocetne matrice A

        finding_min = np.array([])
        indexes = np.array([])
        print(f"Bazne promenljive su : {base}")
        non_base = [(index + 1) for (index, i) in enumerate(base) if i == 0]
        for j in range(1, c.shape[0] + 1):
            # moguce nove bazne promenljive
            if j in non_base:
                temp = w.dot(A[:, j - 1]) - c[j - 1]
                finding_min = np.append(finding_min, temp)
                indexes = np.append(indexes, j)
        
        # od ovih vrednosti treba pronaci najmanji i 
        # takodje treba proveriti da li je on > 0, jer ako je najmanji > 0 onda je algoritam gotov
        min_index = np.argmin(finding_min)
        if max: 
            if finding_min[min_index] > 0:
                print("Najmanji element je > 0, algoritam je gotov")
                print(f"Konacna vrednost max z = {table[0, -2]}")
                print("Maksimalne vrednosti promenljivih:")
                s = 0
                for (index, i) in enumerate(base):
                    if i == 1:
                        s += 1
                        print(f"\tPromenljiva x_{index + 1} = {table[s, -2]}")
                    else: 
                        print(f"\tPromenljiva x_{index + 1} = 0")
                return 
        else:
            if finding_min[min_index] < 0:
                print("Najmanji element je < 0, algoritam je gotov")
                print(f"Konacna vrednost miz z = {-table[0, -2]}")
                print("Maksimalne vrednosti promenljivih:")
                s = 0
                for (index, i) in enumerate(base):
                    if i == 1:
                        s += 1
                        print(f"\tPromenljiva x_{index + 1} = {-table[s, -2]}")
                    else: 
                        print(f"\tPromenljiva x_{index + 1} = 0")
                return 
        
        # ako nije > 0, onda ta promenljiva treba da zameni baznu promenljivu
        # koja je to promenljiva: 
        variable_to_replace = int(indexes[min_index])
        print(f"x_{variable_to_replace} je nova bazna promenljiva")

        # promenljiva koja se zamenjuje dobija se nakon sto se elementi slobodne kolone tackasto podele sa elementima pivot kolone, sto znaci treba postaviti pivot kolonu sledecu
        # pivot kolona racuna se kao: a_j @ B_inv
        a_j = A[:, variable_to_replace - 1]
        print(f"Elementi od kojih nastaje pivot kolona/ koeficijetni koji stoje uz promenljivu x_{variable_to_replace} : {a_j} ")
        pivot_column = B_inv @ a_j
        table[1:, -1] = pivot_column

        # takodje je neophodno staviti slobodan clan u pivot kolonu
        table[0, -1] = finding_min[min_index]

        # sledece treba pronaci promenljivu koju treba zameniti sa novom baznom
        # dot division sa pivot kolonom
        min_dot_div = np.array([])
        for i in range(m): 
            temp = table[i + 1, -2] / table[i + 1, -1]
            min_dot_div = np.append(min_dot_div, temp)

        # najmanji od ovih elemenata, njegova pozicija odredjuje koja ce promenljiva biti zamenjena
        index_of_variable_to_be_replaced = np.argmin(min_dot_div) + 1
        variable_to_be_replaced = int(table[index_of_variable_to_be_replaced, 0])
        print(f"x_{variable_to_be_replaced} je bazna promenljiva koja se menja")

        # neophodno je sada u nultu kolonu na tom indexu staviti naziv promenljive koja je menja
        # kao i u base-vektoru treba da se promeni vrednost na ta dva index-a sa 0 na 1 i sa 1 na 0
        table[index_of_variable_to_be_replaced, 0] = variable_to_replace
        base[variable_to_be_replaced - 1] = 0
        base[variable_to_replace - 1] = 1

        # sada se menjaju redovi u tabeli
        # prvo menjam koja je promenljiva postala bazna, a koja ne bazna
        for i in range(table.shape[0]):
            if i != index_of_variable_to_be_replaced:
                for k in range(table.shape[1] - 2):
                    table[i, k + 1] = table[i, k + 1] - table[index_of_variable_to_be_replaced, k + 1] * table[i, -1] / table[index_of_variable_to_be_replaced, -1]

        table[index_of_variable_to_be_replaced, 1:-1] = table[index_of_variable_to_be_replaced, 1:-1] / table[index_of_variable_to_be_replaced, -1]
        B_inv = table[1:, 1:-2]
    

c = np.array([6, 14, 13])
A = np.array([[0.5, 2, 1], [1, 2, 4]])
b = np.array([[24], [60]])
simplex(c = c, A = A, b = b)

# primer sa min, gde kao bazne promenljive ne uzimamo dodatne vec neke druge (ne x_3, x_4, x_5 nego x_1, x_2, x_4)
# c = np.array([24, 60])
# A = np.array([[-0.5, -1], [-2, -2], [-1, -4]])
# b = np.array([[6], [14], [13]])
# simplex(c = c, A = A, b = b, max = False)

c = np.array([2, 1.5])
A = np.array([[6, 3], [75, 100]])
b = np.array([[1200], [25000]])
simplex(c = c, A = A, b = b, max = True)