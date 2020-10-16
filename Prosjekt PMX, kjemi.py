
####################################################################################################
# Innhold
####################################################################################################

# 1. Konstanter
# 2. Initialisering
# 3. Støkiometri
# 4. Molekyltolkning
# 5. Tabell molmasse
# 6. Tabell reaksjonsligning
# 7. Løse ligninger
# 8. Syrer/baser
# 9. Behandling av input
# 10. Main

####################################################################################################

import string
import sys

####################################################################################################
# 1. KONSTANTER
####################################################################################################

MOLMASSE = 0
MOL = 1
MASSE = 2
KONSENTRASJON = 3
VOLUM = 4
UKJENT = 5


file_data = "pt-data1.csv"
with open(file_data, "r") as f:
    data = f.readlines()

reaklign_string = "Reaksjonsligning"
molmasse_string = "Molmasse"
masser_string = "Masser"
mol_string = "Mol"
konsentrasjoner_string = "Konsentrasjoner"
molekyl_string = "Molekyl"
Ka_str = "Ka"
Kb_str = "Kb"

g_mol = "g/mol"
g = "gram"
mol_l = "mol/L"

pluss = "+"
pluss = f"{pluss:^13}"

pil = "---->"
pil = f"{pil:^17}"

vertikal = "|"
mellomrom = " "
bindestrek = "-"

extra_length = 16


####################################################################################################
# 2. INITIALISERING
####################################################################################################


class Element:
    """klasse for et grunnstoff"""

    def __init__(self, atomnummer, forkortelse, navn, molmasse):
        self.atomnummer = atomnummer
        self.forkortelse = forkortelse
        self.navn = navn
        self.molmasse = molmasse


class Molekyl:
    """klasse for et Molekyl"""

    def __init__(self, molekylformel, molmasse, koeffisient=1, masse=None, konsentrasjon=None, mol=None):
        self.molekylformel = molekylformel
        self.molmasse = molmasse
        self.koeffisient = koeffisient
        self.masse = masse
        self.konsentrasjon = konsentrasjon
        self.mol = mol


class Syre:
    """klasse for Syre"""

    def __init__(self, molekylformel, molekylformel_base, Ka_num, Ka_string, ladning_syre, ladning_base, molmasse_syre, molmasse_base):
        self.molekylformel = molekylformel
        self.molekylformel_base = molekylformel_base
        self.Ka_num = Ka_num
        self.Ka_string = Ka_string
        self.ladning_syre = ladning_syre
        self.ladning_base = ladning_base
        self.molmasse_syre = molmasse_syre
        self.molmasse_base = molmasse_base


def create_element(line):
    """Lager en et objekt til et gitt grunnstoff ut ifra linje i pt-data1.csv"""
    datapoints = line.split(",")
    atomnummer, forkortelse, navn, molmasse = datapoints[:4]

    try:
        molmasse = float(molmasse[:-3])
    except ValueError:
        molmasse = float(molmasse[2:-1])
    return Element(int(atomnummer.strip()), forkortelse.strip(), navn.strip(), molmasse)


def create_molekyler(molekylformler):
    """Lager liste med klasser av molekyler basert på en liste med molekylformler"""
    klasse_molekyler = []
    molmasser, koeffisienter = molmasse_koeff_molekyler(molekylformler)

    for i, molekyl in enumerate(molekylformler):
        if molekyl[0] in string.digits:
            molekylformler[i] = molekyl[1:]

        molekyl = Molekyl(molekylformler[i], molmasser[i], koeffisienter[i])
        klasse_molekyler.append(molekyl)

    return klasse_molekyler


def create_syrer(f_syrer, f_baser, Kaer_nums, Kaer_strings, ladning_syrer, ladning_baser):
    """Lager liste med syreobjekter ut ifra lister med data"""
    syrer = []
    molmasse_syrer = molmasse_koeff_molekyler(f_syrer)[0]
    molmasse_baser = molmasse_koeff_molekyler(f_baser)[0]
    for i in range(len(f_syrer)):
        syre = Syre(f_syrer[i], f_baser[i], Kaer_nums[i],
                    Kaer_strings[i], ladning_syrer[i], ladning_baser[i], molmasse_syrer[i], molmasse_baser[i])
        syrer.append(syre)

    return syrer


def finn_grunnstoff(forkortelse, data):
    """Finner grunnstoff med en gitt forkortelse fra data og gir et objekt"""
    for line in data:
        element = create_element(line)
        if element.forkortelse == forkortelse:
            return element


def check_ferdig(tekst):
    """Avslutter programmet dersom brukeren skriver 'ferdig', mer i bruk ved et program med inputfunksjoner."""
    if tekst == "ferdig":
        print()
        print('Fullført')
        sys.exit()

####################################################################################################
# 3. STØKIOMENTRI
####################################################################################################


def konsentrasjon(V, n):
    c = n / V
    return c


def molmengde1(V, c):
    n = V * c
    return n


def molmengde2(m, molmasse):
    n = m / molmasse
    return n


def masse(n, molmasse):
    m = n * molmasse
    return m


####################################################################################################
# 4. Molekyltolkning
####################################################################################################


def tolke_reaklign(reaklign):
    """Gir lister over reaktanter, produkter og alle molekylene"""
    reaklign = reaklign.replace(" ", "")
    sider = reaklign.split("=")
    reaktanter = sider[0].split("+")
    produkter = sider[1].split("+")
    molekylformler = reaktanter + produkter

    return molekylformler, reaktanter, produkter


def fullfør(molekyl, grunnstoff, liste, i):
    """Sjekker om man er ved slutten på et grunnstoff og legger det til i lista hvis ja"""
    try:
        if molekyl[i + 1] in string.ascii_uppercase or molekyl[i + 1] in "()":
            liste.append(grunnstoff)
    except IndexError:
        liste.append(grunnstoff)


def elements_in_molecule(molekyl):
    """Legger forkortelsene til grunnstoffene i et gitt molekyl til i en
    liste så mange ganger som grunnstoffet er i molekylet og finner eventuell koeffisient"""
    parentes = False
    koeffisient = 1
    grunnstoff_list = []
    list_i_parentes = []
    grunnstoff = ""
    for i in range(len(molekyl)):
        if molekyl[i] in string.ascii_uppercase:
            grunnstoff = molekyl[i]
            fullfør(molekyl, grunnstoff, list_i_parentes, i) if parentes else fullfør(molekyl, grunnstoff, grunnstoff_list, i)

        if molekyl[i] in string.ascii_lowercase:
            grunnstoff += molekyl[i]
            fullfør(molekyl, grunnstoff, list_i_parentes, i) if parentes else fullfør(molekyl, grunnstoff, grunnstoff_list, i)

        if molekyl[i] in string.digits:
            if i == 0:
                koeffisient = molekyl[i]
            elif molekyl[i - 1] == ")":
                pass
            else:
                for a in range(int(molekyl[i])):
                    if parentes:
                        list_i_parentes.append(grunnstoff)
                    else:
                        grunnstoff_list.append(grunnstoff)

        if molekyl[i] in "(":
            parentes = True

        if molekyl[i] in ")":
            x = 1
            parentes = False
            try:
                if molekyl[i + 1] in string.digits:
                    x = int(molekyl[i + 1])
            except IndexError:
                pass
            grunnstoff_list += list_i_parentes * x
    return grunnstoff_list, koeffisient


def molmasse_koeff_molekyler(molekylformler):
    """Gir en liste av molmassene og koeffisientene som tilsvarer molekylformlene i inputlisten"""
    molmasser = []
    koeffisienter = []
    for molekyl in molekylformler:
        grunnstoff_list, koeffisient = elements_in_molecule(molekyl)
        molmasse = 0
        for element in grunnstoff_list:
            grunnstoff = finn_grunnstoff(element, data)
            molmasse += grunnstoff.molmasse
        molmasser.append(molmasse)
        koeffisienter.append(koeffisient)
    return molmasser, koeffisienter


####################################################################################################
# 5. Tabell molmasse
####################################################################################################


def tabellstart():
    """Toppen av tabellen"""
    print("  Molekyl       Molmasse")
    print("--------------------------")


def tabellslutt():
    """Bunnen av tabellen"""
    print("--------------------------")


def tabellelement(molekyl, molmasse):
    """Element i tabellen"""
    print(f"|  {molekyl:14}{molmasse:<8}|")


def tabell_molmasse(molekyler):
    """Lager en formatert tabell av en liste av molmasser og en liste av molekyler som er like lange"""
    molmasser = []
    for molekyl in molekyler:
        molmasser.append(molekyl.molmasse)
    molmasser = list(dict.fromkeys(molmasser))  # Fjerner duplikater

    molekylformler = []
    for molekyl in molekyler:
        molekylformler.append(molekyl.molekylformel)
    molekylformler = list(dict.fromkeys(molekylformler))

    tabellstart()
    for i in range(len(molekylformler)):
        tabellelement(molekylformler[i],
                      round(molmasser[i], 3))
    tabellslutt()


####################################################################################################
# 6. Tabell reaksjonsligning
####################################################################################################


def reaklign_side(side, slutt, reaktanter, produkter):
    """Tar en liste med reaktanter eller produkter og formaterer dem til en string som brukes i tabellen"""
    print_statements = []
    slutter = []
    for molekyl in side:
        if molekyl.koeffisient == 1:
            if molekyl == side[-1]:
                print_statements.append(f"{molekyl.molekylformel}")
                slutter.append(slutt)
            else:
                print_statements.append(f"{molekyl.molekylformel}")
                slutter.append(pluss)
        else:
            if molekyl == side[-1]:
                print_statements.append(
                    f"{molekyl.koeffisient}{molekyl.molekylformel}")
                slutter.append(slutt)
            else:
                print_statements.append(
                    f"{molekyl.koeffisient}{molekyl.molekylformel}")
                slutter.append(pluss)
    final_statement = ""
    for i, statement in enumerate(print_statements):
        if side[i] == reaktanter[-1]:
            print_statements[i] = statement + pil
        else:
            print_statements[i] = statement + slutter[i]
        final_statement += print_statements[i]

    return final_statement, print_statements


def reaklign_final(reaktanter, produkter):
    """Legger sammen de formaterte produktene og reaktantene til en hel, formatert reaksjonsligning"""
    reaktanter_final, reaktanter_statements = reaklign_side(
        reaktanter, "", reaktanter, produkter)
    produkter_final, produkter_statements = reaklign_side(
        produkter, "", reaktanter, produkter)
    reaklign_final = reaktanter_final + produkter_final
    reaklign_final = f"{reaklign_final:^{len(reaklign_final) + extra_length}}"
    return reaklign_final, reaktanter_statements, produkter_statements


def statements_tall(molekyler, reaklign_statements, typ):
    """Printer molmassene i reaksjonsligningstabellen"""
    statement = ""
    for i, molekyl in enumerate(molekyler):
        typer = [molekyl.molmasse, molekyl.mol,
                 molekyl.masse, molekyl.konsentrasjon]
        try:
            tall_string = str(round(typer[typ], 3))
        except TypeError:
            tall_string = "-"
        if i == len(molekyler) - 1:
            statement += f"{tall_string:{len(reaklign_statements[i]) + int(extra_length / 2)}}"
        else:
            statement += f"{tall_string:{len(reaklign_statements[i])}}"

    return statement


def statements_benevning(reaklign_statements, benevning):
    """Printer benevningen direkte under tallet"""
    statement = ""
    for i, molekyl in enumerate(reaklign_statements):
        if i == len(reaklign_statements) - 1:
            statement += f"{benevning:{len(reaklign_statements[i]) + int(extra_length / 2)}}"
            # Legger til ekstra lengde dersom det er det siste statementet i reaksjonsligningen
        else:
            statement += f"{benevning:{len(reaklign_statements[i])}}"

    return statement


def tom_linje(reaklign_output):
    """Printer en tom linje med lengden til reaksjonsligningen og vertikale streker på begge sider"""
    print(f"{vertikal:{len(reaklign_output) + 1}}" + vertikal)


def tabellstrek(reaklign_output, overskrift):
    """Printer en linje med bindestreker og en overskrift i midten"""
    print(f"{overskrift:-^{len(reaklign_output) + 2}}")


def oppdeling(reaklign_output, overskrift):
    """Printer oppdelingene med overskrift"""
    tom_linje(reaklign_output)
    tabellstrek(reaklign_output, overskrift)


def print_benevning(reaklign_statements, typ):
    """Printer benevninger"""
    benevninger = [g_mol, None, g, mol_l]

    print(vertikal + f"{mellomrom:{extra_length / 2}}", end="")
    print(statements_benevning(
        reaklign_statements, benevninger[typ]) + vertikal)


def reaklign_segment(reaklign_output, overskrift):
    """Printer segmentet med reaksjonsligningen"""
    tabellstrek(reaklign_output, reaklign_string)
    tom_linje(reaklign_output)

    print(vertikal + reaklign_output + vertikal)


def tabell2_segment(reaklign_statements, reaklign_output, molekyler, typ):
    """Printer et segment av reaksjonsligninstabellen med en spesiell type info"""
    overskrifter = [molmasse_string, mol_string,
                    masser_string, konsentrasjoner_string]

    oppdeling(reaklign_output, overskrifter[typ])
    tom_linje(reaklign_output)

    print(vertikal + f"{mellomrom:{extra_length / 2}}", end="")
    print(statements_tall(molekyler, reaklign_statements, typ) + vertikal)

    if typ != MOL:
        print_benevning(reaklign_statements, typ)


def reaklign_tabell(reaktanter, produkter, modus):
    """Printer formatert reaksjonsligning og legger reaktantene med slutt og produktene med slutt til i lister"""

    molekyler = reaktanter + produkter
    reaklign_output, reaktanter_statements, produkter_statements = reaklign_final(
        reaktanter, produkter)
    reaklign_statements = reaktanter_statements + produkter_statements

    reaklign_segment(reaklign_output, reaklign_string)

    tabell2_segment(reaklign_statements, reaklign_output, molekyler, MOLMASSE)
    tabell2_segment(reaklign_statements, reaklign_output, molekyler, MASSE)
    tabell2_segment(reaklign_statements, reaklign_output, molekyler, MOL)

    # Printer et segment med konsentrasjoner dersom det er oppgitt for minst ett av molekylene
    if KONSENTRASJON in modus:
        tabell2_segment(reaklign_statements, reaklign_output,
                        molekyler, KONSENTRASJON)

    oppdeling(reaklign_output, bindestrek)


####################################################################################################
# 7. Løse ligninger
####################################################################################################


def masse_addert(molekyler):
    """Adderer massene til en liste med klasser av molekyler"""
    masse_addert = 0
    for molekyl in molekyler:
        try:
            masse_addert += molekyl.masse
        except TypeError:
            pass
    return masse_addert


def differanse(reaktanter, produkter):
    """Regner ut differanse mellom en liste med klasser av reaktanter og produkter"""
    masse_reaktanter = masse_addert(reaktanter)
    masse_produkter = masse_addert(produkter)

    differanse = round(masse_reaktanter - masse_produkter, 3)
    return differanse


def ukjente_molekyler(molekyler):
    """Finner ukjente molekyler"""
    ukjente_molekyler = []
    for molekyl in molekyler:
        if molekyl.masse is None:
            ukjente_molekyler.append(molekyl)
    return ukjente_molekyler


def finn_x_masse(reaktanter, produkter):
    """Finner ukjent masse i to lister av molekylklasser der et eller flere av molekylene har molekyl.masse = None"""
    ukjente_reaktanter = ukjente_molekyler(reaktanter)
    ukjente_produkter = ukjente_molekyler(produkter)

    koeffisient_x = abs(len(ukjente_reaktanter) - len(ukjente_produkter))
    massedifferanse = abs(differanse(reaktanter, produkter))
    try:
        x = massedifferanse / koeffisient_x
    except ZeroDivisionError:
        x = None

    return x


####################################################################################################
# 8. Syrer/baser
####################################################################################################


def f_molekylformler(molekylformler):
    """Tar en liste med molekylformler fra Data Ka.txt og gir en formatert versjon og ladningene"""
    f_molekylformler = []
    ladninger = []
    for molekyl in molekylformler:
        molekyl = molekyl.replace("\n", "")
        molekyl = molekyl.strip()
        if molekyl[-1] in "+-":
            if molekyl[-2] in string.digits:
                molekyl = molekyl.replace(" ", "")
                ladninger.append(molekyl[-2:])
                molekyl = molekyl[:-2]

            else:
                molekyl = molekyl.replace(" ", "")
                ladninger.append(molekyl[-1:])
                molekyl = molekyl[:-1]

        else:
            ladninger.append(None)

        f_molekylformler.append(molekyl)

    return f_molekylformler, ladninger


def read_Ka_data(Ka_data):
    """Leser dataen i Data Ka.txt"""
    org_syrer = []
    org_baser = []
    Kaer_org_strings = []

    for i, line in enumerate(Ka_data):
        if i % 10 == 0:
            chunk = Ka_data[i:i + 10]
            Kaer_org_strings.append(chunk[0])
            org_syrer.append(chunk[4])
            org_baser.append(chunk[6])

    return org_syrer, org_baser, Kaer_org_strings


def tolke_Kaer(Kaer_org_strings):
    """Tar en liste med strings fra Data Ka.txt og gir brukbare tall for Ka"""
    Kaer_nums = []
    Kaer_strings = []

    for Ka_string in Kaer_org_strings:
        Kaer_num = Ka_string.replace(" ", "")
        Kaer_num = Kaer_num.replace("\n", "")
        try:
            tall, eksponentfaktor = Kaer_num.split("*")
            eksponent = float(str(eksponentfaktor)[2:])
            eksponentfaktor = 10 ** eksponent
            Ka_num = float(tall) * eksponentfaktor
            Kaer_nums.append(Ka_num)
            Kaer_strings.append(f"{Ka_num:.3e}")
        except ValueError:
            Kaer_nums.append(None)
            Kaer_strings.append("-")

    return Kaer_nums, Kaer_strings


def syrer_in_reak(syrer, molekyler):
    """Gir syrene i reaksjonsligingen som liste dersom det er noen der"""
    syrer_in_reaklign = []
    mulig_syrereak = False
    for molekyl in molekyler:
        for syre in syrer:
            if molekyl.molekylformel == syre.molekylformel:
                if molekyl.molekylformel not in ("H2O", "H3O"):
                    mulig_syrereak = True
                    syrer_in_reaklign.append(syre)

    return syrer_in_reaklign, mulig_syrereak


def mulige_syrer():
    """Gir en liste over alle syrene i Data Ka.txt som objekter"""
    with open("Data Ka.txt", "r") as f:
        Ka_data = f.readlines()

    org_syrer, org_baser, Kaer_org_strings = read_Ka_data(Ka_data)

    Kaer_nums, Kaer_strings = tolke_Kaer(Kaer_org_strings)

    f_syrer, ladning_syrer = f_molekylformler(org_syrer)
    f_baser, ladning_baser = f_molekylformler(org_baser)

    syrer = create_syrer(f_syrer, f_baser, Kaer_nums,
                         Kaer_strings, ladning_syrer, ladning_baser)

    return syrer


def print_syresegment(syre, Kb, padding, k, len_tabell):
    print(f"{vertikal:{padding}}", end="")
    print(
        f"{syre.molekylformel:{k}}{round(syre.molmasse_syre, 3):<{k}}{syre.Ka_string:{k}}{bindestrek:{padding + 1}}", end="")
    print(f"{vertikal:}")
    print(f"{vertikal:{padding}}", end="")
    print(
        f"{syre.molekylformel_base:{k}}{round(syre.molmasse_base, 3):<{k}}{bindestrek:{k}}{Kb:{padding + 1}}", end="")
    print(f"{vertikal:}")


def syretabell(syrer, molekyler):
    """Tar en liste med syrer og printer tabell med molekylformel, molmasse, Ka og Kb til korresponderende base"""
    Kw = 10e-14
    k = 16
    padding = 13
    print(f"{mellomrom:{padding}}{molekyl_string:{k}}{molmasse_string:{k}}{Ka_str:{k}}{Kb_str:{padding}}")
    # 4 for starten, 3k for alt til og med Ka, 2 for Kb og 4 for symmetri
    len_tabell = padding + k * 3 + 2 + padding
    print(bindestrek * len_tabell)
    for syre in syrer:
        try:
            Kb = f"{Kw / syre.Ka_num:.3e}"
        except TypeError:
            Kb = bindestrek
        print_syresegment(syre, Kb, padding, k, len_tabell)
    print(bindestrek * len_tabell)


####################################################################################################
# 9. Behandling av inpput
####################################################################################################


def gram(svar):
    if "kg" in svar:
        return float(svar.strip("kg")) * 1000
    else:
        return float(svar.strip("g"))


def liter_til_gram(svar, c, molekyler, i):
    """Regner ut massen til et molekyl ved å vite molekylet, volum og konsentrasjon"""
    if "ml" in svar:
        svar = V = float(svar.strip("ml")) / 1000
    else:
        svar = V = float(svar.strip("l"))
    molmengde = molmengde1(V, c)
    return round(float(masse(molmengde, molekyler[i].molmasse)), 3)


def konsentrasjon_til_gram(svar, V, molekyler, i):
    """Omgjør mol/l til gram ved å vite volumet"""
    c = float(svar.strip("mol/l"))
    molmengde = molmengde1(V, c)
    return round(masse(molmengde, molekyler[i].molmasse), 3)


def tolke_input(inputs, molekyler):
    """Tolker input av masse, liter, ukjent eller konsentrasjon til de respektive atomene"""
    masser = []
    volum = None
    modus = []
    for i, svar in enumerate(inputs):
        svar = svar.replace(" ", "").lower()
        try:
            if svar in "x":
                modus.append(UKJENT)
                masser.append(None)

            if "g" in svar:
                modus.append(MASSE)
                masser.append(gram(svar))

            if svar[-1] not in string.ascii_lowercase:
                modus.append(MOL)
                masser.append(float(svar) * molekyler[i].molmasse)

            if "l" in svar and "mol/l" not in svar:
                modus.append(KONSENTRASJON)
                c = 0.5  # uferdig for øyeblikket, i en fullstendig versjon kan man la brukeren bestemme konsentrasjonen
                masser.append(liter_til_gram(svar, c, molekyler, i))
                molekyler[i].konsentrasjon = c

            if "mol/l" in svar:
                modus.append(KONSENTRASJON)
                V = 2  # uferdig for øyeblikket, i en fullstendig versjon kan man la brukeren bestemme konsentrasjonen
                masser.append(konsentrasjon_til_gram(svar, V, molekyler, i))
                c = float(svar.strip("mol/l"))
                molekyler[i].konsentrasjon = c
        except IndexError:
            pass
    return masser, volum, modus


####################################################################################################
# 10. MAIN
####################################################################################################

#
# Programmet kan:
#
#     - Regne ut ukjente masser i en fullstending reaksjonsligning og vise diverse informasjon i en formatert tabell.
#     - Vise en tabell med molmasser ut ifra en liste med molekyler.
#     - Vise en tabell med Ka og Kb dersom den finner ut at et eller flere av molekylene, enten i reaksjonsligningen eller i
#       listen med molekyler, er en syre.
#
# Endre områdene som er formatert på måten under for å bruke programmet.
##########################
# Instruksjoner
# variabel = liste/string
##########################



def main():

    ####################################################################################################################
    # Endre dette for å aktivere de ulike delene av programmet, den siste bør alltid være "ferdig" for å unngå at programmet crasher:
    initinputs = ["ukjent", "molmasser", "ferdig"]  # Muligheter: "ukjent", "molmasser" og "ferdig". Dette er laget som eksempel på hvordan man kan gjøre det med inputs senere
    ####################################################################################################################

    for n in range(100):  # Kan ha en while loop her, men er litt ork hvis det er en bug og den holder på for alltid

        mulig_syrereak = False
        ukjent_masse = False
        molmasser_liste = False

        initinput = initinputs[n].lower()

        if initinput == 'ukjent':
            ukjent_masse = True

        if initinput == 'molmasser':
            molmasser_liste = True

        check_ferdig(initinput)
        if ukjent_masse:
            #######################################################################################
            # Endre dette dersom du ønsker å finne ukjent masse for forskjellige molekyler:
            reaklign = "HNO2 + 4NH3 = Cu(NH3)4 + Ba(OH)2"
            #######################################################################################

            molekylformler, reaktanter, produkter = tolke_reaklign(reaklign)
            molmasser, koeffisienter = molmasse_koeff_molekyler(
                molekylformler)  # Finner molmasser og koeffisienter)
            molekyler = create_molekyler(molekylformler)

            ###########################################################################################################################################
            # Endre dette dersom du ønsker å finne ukjent masse for forskjellige molekyler
            # Benevninger: ingenting (mol), g/kg (gram/kilogram), l/ml (volum), mol/l (konsentrasjon)
            inputs = ["x", "x", "4 mol/l", "2"]  # Default konsentrasjon er 0.5 og default volum er 2, dette kan endres på linje 689 og 696
            # Pass på at det er like mange eller flere inputs som molekyler for å unngå at programmet krasjer
            ###########################################################################################################################################

            masser, volum, modus = tolke_input(inputs, molekyler)

            for i, molekyl in enumerate(molekyler):
                molekyl.masse = masser[i]

            reaktanter = molekyler[:len(reaktanter)]
            produkter = molekyler[len(reaktanter):]

            x = finn_x_masse(reaktanter, produkter)
            ukjente_molekyler1 = ukjente_molekyler(molekyler)

            if x is not None:
                for molekyl in (ukjente_molekyler1):
                    molekyl.masse = round(x, 3)

            for molekyl in molekyler:
                try:
                    molekyl.mol = molekyl.masse / molekyl.molmasse
                except TypeError:
                    pass

            reaklign_tabell(reaktanter, produkter, modus)
            print()
            print()
            print()

        # Brukeren ønsker ikke å bruke reaksjonsliginger, men bare vil ha en liste over molmasser
        if molmasser_liste:
            try:

                #######################################################################################################
                # Endre dette dersom du ønsker molmasser for forskjellige molekyler, separer molekylene med et komma:
                molecules = "BaSO4, Ca, SO4, (COOH)2, HCl, HNO3"
                #######################################################################################################

                tolkbart_input = molecules.replace(" ", "")
                molekylformler1 = tolkbart_input.split(",")
                molekyler = create_molekyler(molekylformler1)
                tabell_molmasse(molekyler)
                print()
                print()
                print()
            except AttributeError:
                pass

        syrer = mulige_syrer()
        syrer_in_reaklign, mulig_syrereak = syrer_in_reak(syrer, molekyler)

        if mulig_syrereak:
            syretabell(syrer_in_reaklign, molekyler)
            print()
            print()
            print()


if __name__ == '__main__':
    main()
