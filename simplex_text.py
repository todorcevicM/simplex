"""
max z = 6x_1 + 14x_2 + 13x_3

0.5x_1 + 2x_2 + 1x_3 <= 24
1x_1 + 2x_2 + 4x_3 <= 60
x_1, x_2, x_3 >= 0

1. ukoliko su neke od jednakosti >= mora da se pomnozi ta jednakost sa -1
2. pretvaramo <= u < tako sto dodajemo nove promenljive: 
    0.5x_1 + 2x_2 + 1x_3 + x_4 = 24
    1x_1 + 2x_2 + 4x_3 + x_5 = 60
3. dobijamo vektor c = [6, 14, 13, 0, 0]
4. dobijamo matricu 
A_1 = 
    [
        [0.5, 2, 1, 1, 0], 
        [1  , 2, 4, 0, 1]
    ]
4.1. dobijamo matricu b = [24, 60]^T = [[24], [60]]

i sad zapravo krecu koraci revidiranog simplex-a
1. odredjivanje baznih promenljivih, baznih promenljivih ima koliko ima ogranicenja: A_1.shape[0] = 2
    dodatno ogranicenje pri biranju baznih promenljivih jeste da matrica napravljena od koeficijenata tih promenljivih mora imati inverznu matricu
    dodatno ogranicenje je takodje da ta inverzna matrica pomnozena sa matricom b treba da dobijemo brojeve >= 0
    najjednostavnije resenja za biranje baznih promenljivih su uglavnom dodatne promenljive x_4 i x_5
    znaci za matricu baznih promenljivih imamo B = [[1, 0], [0, 1]]
    u tom slucaju inverzna matrica B^-1 = x_B ce biti ista kao i B: B^-1 = [[1, 0], [0, 1]]
    B^-1 * b = [[1, 0], [0, 1]] * [[24], [60]] = [[24], [60]]

2. pravimo tabelu, u prvu kolonu stavljamo bazne promenljive (x_4 i x_5), u slucaju algoritma mogu biti (indeksi + 1) (indeks krece od 0)
    u sledecih n kolona stavljamo B^-1 matricu, 
    prvi red ostaje prazan i sada se popunjava sa w (omega) w = c_B * B^-1, gde je c_B vektor koji sadrzi koeficijente baznih promenljivih onome sto zelimo da maksimizujemo, u ovom slucaju c_B = [0, 0]
    u sledecu slobodnu kolonu (bez da stavljamo u prvi red) stavljamo b^- = B^-1 * b, b^- = [[1, 0], [0, 1]] * [[24], [60]] = [[24], [60]]
    u prvi red, u istoj toj koloni ide vrednost koja se dobija kao: c_B * b^- = [0, 0] * [[24], [60]] = 0

    sledeca kolona jeste pivot kolona, i odredjuje se na sledeci nacin:
    za sve ne bazne promenljive (x_1, x_2, x_3) racunamo: u ovom slucaju bice korisceni j = (indeks + 1) 
        za j: w * a_j - c_j
        gde je a_j transponovani vektor koji u sebi sadrzi kolonu sa koeficijentima promenljivih, u toj koloni (zapravo ce biti j - 1) 
            znaci prolazi kroz sva ogranicenja i uzima koeficijente promenljivih koje se nalaze u koloni j - 1
        i c_j je jedan element vektora c, koji sadrzi koeficijent koji se nalazi u koloni j - 1 (odnosno na mestu j - 1)
        j = 1: w * a_1 - c_1 = [0, 0] * [[0.5], [1]] - [6] = [-6]
        j = 2: w * a_2 - c_2 = [0, 0] * [[2], [2]] - [14] = [-14]
        j = 3: w * a_3 - c_3 = [0, 0] * [[1], [4]] - [13] = [-13]
    sada mozemo da odaberemo novu baznu promenljivu, nju dobijamo tako sto cemo odabrati najnegativniji element, u ovom slucaju za j = 2, jer je -14 najnegativniji element
    to znaci da ce nam u sledecoj iteraciji bazna promenljiva biti x_2

    sledece popunjavamo poseldnju kolonu tabele
        u prvi red te kolone upisujemo tu najnegativniju vrednost -14, 
        u ostatak kolone stavljamo a_j * B^-1, u ovom slucaju a_2 * B^-1 = [[2], [2]] * [[1, 0], [0, 1]] = [[2], [2]]
        vidimo da su nam obe ove vrednosti pozitivne i to znaci da algoritam jos uvek nije gotov
    dobili smo pivot kolonu
    pivot kolona ce se nalaziti na poziciji: n + 2
    sledeci trazimo pivot vrstu koja se dobija tako sto se vrednosti slobodne kolone podele sa odgovarajucim vrednostima pivot kolone
    tackasto deljenje matrica: [[24], [60]] / [[2], [2]] = [[12], [30]]
    kao pivot kolonu biramo onu koja ima najmanju vrednost od dobijenih vrednosti koje smo dobili deljenjem [[24], [60]] / [[2], [2]] i u ovom slucaju je to 12, sto znaci prvi red p_v
    
    pivot element dobija se u preseku pivot kolone i pivot vrste i u matrici koju cu praviti se nalazi na poziciji [j][n + 2] = 2

        |   0   0   |   0   |   -14 |
    --------------------------------|
    x_4 |   1   0   |   24  |   2   |
    x_5 |   0   1   |   60  |   2   |   

3. menjamo taj j-ti red sa novim redom koji ce se racunati
    u nultu kolonu dolazi oznaka elementa koji se dodaje, x_2
    dok se ostali elemetni tog reda dobijaju tako sto se prethodne vrednosti dele sa pivot elementom
    [1, 0] / 2 = [0.5, 0]
    i sledeca, n + 1 kolona se racuna na isti nacin i to ce ovde biti [24] / 2 = [12]

    ostali redovi se racunaju na drugaciji nacin, i to tako: 
        (taj element) - ((element iz pivot vrste koji se nalazi u toj koloni) * (element iz pivot kolone koji se nalazi u tom redu)) / (pivot element)
        i to je:
            0 - (1 * 2) / 2 = -1

        for (i in 1 to n):
            if (i != j):
                for (k in 1 to n):
                    a[i][k] = a[i][k] - (a_[j][k] * a[k][n + 2]) / a[j][n + 2]
        to se isto radi i za n + 1 kolonu

    na kraju dobijamo: 

        |   7    0   |   168 |   -6    |
    -----------------------------------|
    x_2 |   0.5  0   |   12  |   0.5   |
    x_5 |   -1   1   |   36  |   3     | 

    kao novo B^-1 matricu: B^-1 = [[0.5, 0], [-1, 1]]
    takodje, c_B postaje: c_B = [14, 0] 

    ovde smo izracunali w ponovo: w = c_B * B^-1 = [14, 0] * [[0.5, 0], [-1, 1]] = [7, 0]

    ponovo racunamo za j = 1, 3, 4
    za j = 1: w * a_1 - c_1 = [7, 0] * [[0.5], [1]] - [6] = -2.5
    za j = 3: w * a_3 - c_3 = [7, 0] * [[1], [4]] - [13] = -6
    za j = 4: w * a_4 - c_4 = [7, 0] * [[1], [0]] - [0] = 7
    opet biramo najnegativniji element, u ovom slucaju za j = 3, jer je -6 najnegativniji element
    samim tim x_3 ce biti nova bazna promenljiva

    sledece popunjavamo poseldnju kolonu tabele
        u prvi red te kolone upisujemo tu najnegativniju vrednost -6,
        u ostatak kolone stavljamo a_j * B^-1, u ovom slucaju a_3 * B^-1 = [[1], [4]] * [[0.5, 0], [-1, 1]] = [[0.5], [3]]
        vidimo da su nam obe ove vrednosti pozitivne i to znaci da algoritam jos uvek nije gotov
    dobili smo pivot kolonu
    pivot kolona ce se nalaziti na poziciji: n + 2
    sledeci trazimo pivot vrstu koja se dobija tako sto se vrednosti slobodne kolone podele sa odgovarajucim vrednostima pivot kolone
    tackasto deljenje matrica, cije su vrednosti uzete iz [n + 1] kolone i [n + 2] kolone: [[12], [36]] / [[0.5], [3]] = [[24], [12]]
    biramo najmanji pozitivan broj i u ovom slucaju je to 12, red u kojem se on nalazi oznavava pivot vrstu

    racunamo ponovo koju cemo staru baznu promenljivu izbaciti iz B^-1 matrice 
    i u nju upisujemo novu baznu promenljivu, u ovom slucaju x_3


        |   5    2      |   240  |  -1.5  |
    --------------------------------------|
    x_2 |   2/3  -1/6   |   6    |  1/6   |
    x_3 |   -1/3  1/3   |   12   |  1/6   | 

    kao novo B^-1 matricu: B^-1 = [[2/3, -1/6], [-1/3, 1/3]]
    takodje, c_B postaje: c_B = [14, 13]

    ovde smo izracunali w ponovo: w = c_B * B^-1 = [14, 13] * [[2/3, -1/6], [-1/3, 1/3]] = [5, 2]

    ponovo racunamo za j = 1, 4, 5
    za j = 1: w * a_1 - c_1 = [5, 2] * [[0.5], [1]] - [6] = -1.5
    za j = 4: w * a_4 - c_4 = [5, 2] * [[1], [0]] - [0] = 5
    za j = 5: w * a_5 - c_5 = [5, 2] * [[0], [1]] - [0] = 2

    sledece popunjavamo poseldnju kolonu tabele
        u prvi red te kolone upisujemo tu najnegativniju vrednost -1.5,
        u ostatak kolone stavljamo a_j * B^-1, u ovom slucaju a_1 * B^-1 = [[0.5], [1]] * [[2/3, -1/6], [-1/3, 1/3]] = [[1/6], [1/6]]
        vidimo da su nam obe ove vrednosti pozitivne i to znaci da algoritam jos uvek nije gotov
    dobili smo pivot kolonu
    pivot kolona ce se nalaziti na poziciji: n + 2
    sledeci trazimo pivot vrstu koja se dobija tako sto se vrednosti slobodne kolone podele sa odgovarajucim vrednostima pivot kolone
    tackasto deljenje matrica, cije su vrednosti uzete iz [n + 1] kolone i [n + 2] kolone: [[6], [12]] / [[1/6], [1/6]] = [[36], [72]]

    racunamo ponovo koju cemo staru baznu promenljivu izbaciti iz B^-1 matrice
    i u nju upisujemo novu baznu promenljivu, u ovom slucaju x_1

        |   11   0.5 |   294   |  -1.5  |
    --------------------------------------|
    x_1 |   4   -1   |   36    |  1/6   |
    x_3 |  -1  0.5   |    6    |  1/6   |   

    kao novo B^-1 matricu: B^-1 = [[4, -1], [-1, 0.5]]
    takodje, c_B postaje: c_B = [6, 13]

    ovde smo izracunali w ponovo: w = c_B * B^-1 = [6, 13] * [[4, -1], [-1, 0.5]] = [11, 0.5]

    ponovo racunamo za j = 2, 4, 5
    za j = 2: w * a_2 - c_2 = [11, 0.5] * [[2], [2]] - [6] = 9
    za j = 4: w * a_4 - c_4 = [11, 0.5] * [[1], [0]] - [0] = 11
    za j = 5: w * a_5 - c_5 = [11, 0.5] * [[0], [1]] - [0] = 0.5

    dobili smo sve pozitivne elemente sto znaci da je ovo kraj algoritma

        algoritam nastavlja sa radom sve dok ne dobijemo sve pozitivne vrednosti kada se trazi nova bazna promenljiva
        ukoliko su sve vrednosti pozitivne, dobijamo resenje

    resenje dobijamo citanjem vrednosti iz slobodne kolone: 
        resenje maksiziranja max z = 294, sa pozicije [0][n + 1]
        vrednosti koeficijenata promenljivih x_1, x_2, x_3, x_4, x_5 su:
            ako se nalaze u tabeli: 
                x_1 = 36 
                x_3 = 6
            ako se ne nalaze u tabeli: 
                x_2 = 0
                x_4 = 0
                x_5 = 0
"""