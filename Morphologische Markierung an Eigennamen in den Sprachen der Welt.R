############### Allgemeiner Teil – Funktion Morph einzeln und für ein ganzes Korpus ############

# Laden der nötigen Pakete
library(qlcMatrix)
library(readr)
library(stringdist)
library(maps)
library(mapproj)

# Berechne für einen Text (NT), wie viele Formen für die Übersetzung eines Namens vorkommen
# param:
# – file_x: (string) Absoluter Dateipfad zu einer Bibel der Zielsprache
# – dateiBibel1 (string) absoluter Dateipfad zur Ausgangsbibel
# – name: (string) der Name, der als Grundlage dient. Default = "petrus"
# – th1: (float) Schwelle zwischen 0 und 1, ab deren Ähnlichkeit Wortformen aufgenommen werden (nur falls method = sim). Default = 0.5
# – th2: (int) Schwelle zwischen 1 und length(wordforms), die die Länge der sortierten assoc-Werte bestimmt
# – method (string) default sim führt zu Stringvergleich per sim.strings, alle anderen werden 
#                   an stringdist weitergereicht
# return int > 0: Anzahl der gefundenen Namenvarianten

morph = function(file_x,
                 dateiBibel1,
                 th1 = 400,
                 th2 = 0.6,  
                 name = "petrus", 
                 method = "sim") {
  
  # Sichergehen, dass der Name in Minuskeln vorliegt
  name = tolower(name)
  
  #Einlesen der Ausgangsbibel und der Zielbibel
  deu = read.text(dateiBibel1)
  x = read.text(file_x)
  
  # Schnittmenge über vorhandene Verse
  globalID = intersect(names(deu),names(x))
  if(length(globalID) == 0) {return(0)}
  
  # Aus den Versstrings werden Wortformenstrings erzeugt, 
  # sodass eine Wort-Vers-Matrix entsteht (lowercase = T)
  WV_deu = splitText(deu, globalSentenceID = globalID, simplify = T)
  WV_x = splitText(x, globalSentenceID = globalID, simplify = T)
  
  # Alle Verse der deutschen Übersetzung, in denen name vorkommt
  # verseMitName = which(WV_deu[name,] > 0)
  # Welche Wörter kommen in den parallelen Versen der Zielsprache absolut am häufigsten vor?
  # sort(rowSums(WV_x[,verseMitName]), decreasing = T)[1:20]
  
  # Besser ist es, die Korrelationen zwischen den Zeilen zu bestimmen
  # dafür zunächst Transformieren in binäre Matrizen
  WV_deu[which(WV_deu > 0)] = 1
  WV_x[which(WV_x > 0)] = 1
  assoc = assocSparse(t(WV_deu), t(WV_x), method = res)
  
  
  # Höchste Korrelationen
  kandidaten = sort(assoc[name,], decreasing = T)[1:th1]
  
  # Als Entsprechung des Namens in der Zielsprache wird
  # vorläufig einfach der mit der höchsten Korrelation gewählt
  compName = names(kandidaten[1])
  
  #Anhand dieser Entsprechung werden nun ähnliche Forme gesucht
  wordforms = rownames(WV_x)
  
  # Wahl der Methode zum Stringvergleich
  if (method =="sim") {
    #Ähnlichkeit zu allen Wortformen
    #liefert Werte zwischen 0 und 1
    sim = sort(sim.strings(compName, wordforms, boundary = T), decreasing = T)
    #Gib die strings der besten (ab einer gewissen Schwelle th aus)
    sim = names(sim[which(sim >= th2)])
  }
  else {
    #Ähnlichkeit zu allen Wortformen
    distMat = stringdistmatrix(compName, rownames(WV_x),
                               method = method,
                               useNames = "strings",
                               weight = c(d = 1, i = 0.5, s = 0.9, t = 1),
                               q = 3)
    sim = sort(distMat[1,])
    
    #Gib die besten strings (ab einer gewissen Schwelle) zurück
    sim = names(sim[which(sim < th2)])
  }
  
  #Schnittmenge aus beiden Listen
  sims = intersect(sim, names(kandidaten))
  
  # Gib die gefundenen Formen auf der Konsole aus
  print(sims)
  
  return(length(sims))
}



####################### Anwenden auf alle Bibeln ###########################

# Führt die Funktion morph(x, y) für alle Bibeln in einem Verzeichnis aus
# die Dateinamen in dem Verzeichnis sollten mit dem ISO-Code der Sprachen anfangen
#param:
# corpus_dir: (string) absoluter Verzeichnispfad des Korpus
# bib1: (string) absoluter Dateipfad der Ausgangsbibel wird an morph weitergereicht
# name: (string) Name, der an morph weitergereicht wird, default = petrus
#return: named numeric mit Sprach-ISO-Codes als Namen und morph-Werten als values
execMorph = function(corpus_dir, bib1, name = "petrus") {

  # Initialisierung
  # Alle Dateipfade des Verzeichnis
  bibles = (list.files(corpus_dir, full.names = T))
  
  # Davon jeweils die ersten drei Buchstaben
  biblecodes = substring(list.files(corpus_dir, full.names = F), 1, 3)
  
  # Alle Anfangsbuchstaben (ISO-Codes) ohne Duplikate
  biblecodesSet = intersect(biblecodes,biblecodes)
  
  # Allokieren und benennen des Ergebnisvektors
  results = rep(0,length(biblecodesSet))
  names(results) = biblecodesSet
  
  # Randbetrachtung aus Indexgründen
  results[1] = morph(bibles[1], bib1, name=name)
  
  if(length(bibles) == 1) {return(results)}
  # Hier wird die Funktion morph auf alle Dateien aus corpus_dir angewendet;
  # falls dabei meherere die gleichen drei Anfangsbuchstaben (ISO-Code) tragen
  # wird der Mittelwert aus diesen gebildet
  k = 2
  z = 0
  for (i in 2:length(bibles)) {
    if(biblecodes[i] == biblecodes[i-1]) {
      results[k-1] = results[k-1] + morph(bibles[i],bib1,name=name)
    }
    else {
      results[k] = morph(bibles[i], bib1, name = name)
      if(i != k + z) {
        results[k-1] = floor(results[k-1]/(i-k+1-z))
        z = i-k
      }
      k = k + 1
    }
  }
  
  return(results)
}





############### Spezifischer Teil – Bearbeitungsschritte und interessante Ausgaben (hardcode!) ################
### Hier werden zusätzliche Interessante Ideen, die nicht von meinem Dateisystem abstrahiert sind,
### aufgelistet. Dazu zählt zum Beispiel, die Optimierung des Algorithmus, die im Dialog
### mit verschiedenen Konsolenausgaben stattgefunden hat. Das Plotten der Karten könnte noch interessant sein
### und ist möglichst nachvollziehbar und allgemein gehalten.



#### Plotten der Karten für Sprachen mit und ohne morphol. Markierung ####

res = execMorph("/home/fabian/zack/Bibelkorpus/", 
                dateiBibel1 = "/home/fabian/zack/Bibelkorpus/deu-x-bible-Genfer2011")
resultsDavid = execMorph("/home/fabian/zack/Bibelkorpus/",
                         dateiBibel1 = "/home/fabian/zack/Bibelkorpus/deu-x-bible-genfer2011.txt",
                         name = "david")

# Auswerten und Einteilen der Ergebnisse
hasMorph = res > 1
hasMorphDavid = resultsDavid > 1

# Zuteilen der ISO-Codes
biblecodesYes = names(res[hasMorph])
biblecodesNo  = names(res[!hasMorph])

# Einlesen der Info Datei, die die Koordinaten zu den Codes enthält
bibleinfo = read_tsv("/home/fabian/zack/ISO639-3.tsv")

# Zuschneiden auf vorhandene ISO-Codes
biblecodesSet = intersect(biblecodes,biblecodes) # Das ist das gleiche wie names(res)

biblecodesYesForMap = intersect(biblecodesYes, bibleinfo$ISO)
biblecodesNoForMap  = intersect(biblecodesNo, bibleinfo$ISO)

# Nicht schön (depricated), aber tut, was es soll
rownames(bibleinfo) = bibleinfo$ISO
# Nämlich extrahieren der Koordinaten
coordYesLON = bibleinfo[biblecodesYesForMap,]$LMP_LON
coordYesLAT = bibleinfo[biblecodesYesForMap,]$LMP_LAT
coordNoLON  = bibleinfo[biblecodesNoForMap,]$LMP_LON
coordNoLAT  = bibleinfo[biblecodesNoForMap,]$LMP_LAT

# Zeichnen der Karte
map()
text(coordYesLON, coordYesLAT, ".", col = "red")
text(coordNoLON, coordNoLAT, ".", col = "red")

text(coordYesLONDavid,coordYesLATDavid, ".", col= "green")
text(coordNoLONDavid, coordNoLATDavid, ".", col = "red")


# Anzahl an Sprachen mit und ohne Morphologie am EN 
lanWithENMorph = length(biblecodesYes)
lanWithoutENMorph = length(biblecodesNo)

#Höchster erreichter Morphwert
max(resultsDavid)




######################### Einlesen einiger Bibeln ###############################
are = "/home/fabian/zack/Bibelkorpus/are-x-bible.txt"
apz = "/home/fabian/zack/Bibelkorpus/apz-x-bible.txt"
nyf = "/home/fabian/zack/Bibelkorpus/nyf-x-bible.txt"
ksw = "/home/fabian/zack/Bibelkorpus/ksw-x-bible.txt"
guw = "/home/fabian/zack/Bibelkorpus/guw-x-bible-newworld.txt"
arb = "/home/fabian/zack/Bibelkorpus/arb-x-bible.txt"
lit = "/home/fabian/zack/Bibelkorpus/lit-x-bible-ecumenical.txt"
lit2= "/home/fabian/zack/Bibelkorpus/lit-x-bible-tikejimozodis.txt"
cme = "/home/fabian/zack/Bibelkorpus/cme-x-bible.txt"
deu = "/home/fabian/zack/Bibelkorpus/deu-x-bible-schlachter2000.txt"
tgl = "/home/fabian/zack/Bibelkorpus/tgl-x-bible-1996.txt"
aso = "/home/fabian/zack/Bibelkorpus/aso-x-bible.txt"
eus = "/home/fabian/zack/Bibelkorpus/eus-x-bible-batua.txt"
aak = "/home/fabian/zack/Bibelkorpus/aak-x-bible.txt"
nhg = "/home/fabian/zack/Bibelkorpus/nhg-x-bible.txt"
ita = "/home/fabian/zack/Bibelkorpus/ita-x-bible-newworld.txt"
zul = "/home/fabian/zack/Bibelkorpus/zul-x-bible-newworld.txt"
zul2 = "/home/fabian/zack/Bibelkorpus/zul-x-bible.txt"
zyp = "/home/fabian/zack/Bibelkorpus/zyp-x-bible.txt"
moh = "/home/fabian/zack/Bibelkorpus/moh-x-bible.txt"
jvn = "/home/fabian/zack/Bibelkorpus/jvn-x-bible.txt"
mri = "/home/fabian/zack/Bibelkorpus/mri-x-bible.txt"
arz = "/home/fabian/zack/Bibelkorpus/arz-x-bible.txt"
nav = "/home/fabian/zack/Bibelkorpus/nav-x-bible.txt"
pan = "/home/fabian/zack/Bibelkorpus/pan-x-bible.txt"
swe = "/home/fabian/zack/Bibelkorpus/swe-x-bible-2000.txt"
fin = "/home/fabian/zack/Bibelkorpus/fin-x-bible-2012.txt"
eng = "/home/fabian/zack/Bibelkorpus/eng-x-bible-montgomery.txt"
pol = "/home/fabian/zack/Bibelkorpus/pol-x-bible-covenant.txt"
zho = "/home/fabian/zack/Bibelkorpus/zho-x-bible-standard.txt"
ewe = "/home/fabian/zack/Bibelkorpus/ewe-x-bible.txt"
cmo = "/home/fabian/zack/Bibelkorpus/cmo-x-bible-khmr.txt"
pck = "/home/fabian/zack/Bibelkorpus/pck-x-bible.txt"
fij = "/home/fabian/zack/Bibelkorpus/fij-x-bible-bau.txt"
arz = "/home/fabian/zack/Bibelkorpus/arz-x-bible.txt"
tur = "/home/fabian/zack/Bibelkorpus/tur-x-bible-2009.txt"
ton = "/home/fabian/zack/Bibelkorpus/ton-x-bible.txt"
qve = "/home/fabian/zack/Bibelkorpus/qve-x-bible.txt"
hau = "/home/fabian/zack/Bibelkorpus/hau-x-bible.txt"
spa = "/home/fabian/zack/Bibelkorpus/spa-x-bible-paratodos.txt"
soy = "/home/fabian/zack/Bibelkorpus/soy-x-bible.txt" # Bible auf Miyobe
hin = "/home/fabian/zack/Bibelkorpus/hin-x-bible-common.txt"
chk = "/home/fabian/zack/Bibelkorpus/chk-x-bible.txt"
cko = "/home/fabian/zack/Bibelkorpus/cko-x-bible.txt"
yor = "/home/fabian/zack/Bibelkorpus/yor-x-bible-2010.txt"
ifk = "/home/fabian/zack/Bibelkorpus/ifk-x-bible.txt"
naq = "/home/fabian/zack/Bibelkorpus/naq-x-bible.txt"
swh = "/home/fabian/zack/Bibelkorpus/swh-x-bible-union1997.txt"
yrk = "/home/fabian/zack/Bibelkorpus/yrk-x-bible.txt"
urd = "/home/fabian/zack/Bibelkorpus/urd-x-bible-devanagari.txt"
urdlat = "/home/fabian/zack/Bibelkorpus/urd-x-bible-latn.txt"
deu2 = "/home/fabian/zack/Bibelkorpus/deu-x-bible-freebible.txt"
deuGenfer = "/home/fabian/zack/Bibelkorpus/deu-x-bible-genfer2011.txt"
deuluther = "/home/fabian/zack/Bibelkorpus/deu-x-bible-luther2017.txt"
deulutherAlt = "/home/fabian/zack/Bibelkorpus/deu-x-bible-luther1545letztehand.txt"

######################## Optimierung von morph(x) #########################
morph(qve, deuGenfer, th1 = 800, th2 = 0.45)

morph(qve, deuGenfer, th1 = 400, th2 = 0.7)

######################### Testen ##########################

# Khoikhoi khoisan-sprache in Namibia
morph(naq, deuGenfer)

# Maori Austronesische Sprache
morph(mri, deuGenfer)

# Türkisch, Turksprache
morph(tur, deuGenfer)

# Polnisch, Indogermanisch
morph(pol, deuGenfer)

# Ewe, Niger-Kongo Sprache
morph(ewe, deuGenfer)

#Telecingo Nahuatl, uto-atztekisch
morph(nhg, deuGenfer)


# Mit anderen Namen

### David ###

# Khoikhoi khoisan-sprache in Namibia
morph(naq, deuGenfer, name = "david")

# Maori, Austronesisch
morph(mri, deuGenfer, name = "david")

# Türkisch
morph(tur, deuGenfer, name = "david")

# Polnisch, Indogermanisch
morph(pol, deuGenfer, name = "david")

# Ewe, Niger-Kongo
morph(ewe, deuGenfer, name = "david")

#Telecingo Nahuatl, uto-atztekisch
morph(nhg, deuGenfer, name = "david")



### Martha ###

# Khoikhoi khoisan-sprache in Namibia
morph(naq, deuGenfer,name = "martha")

# Maori, Austronesisch
morph(mri,deuGenfer, name = "martha")

# Türkisch
morph(tur, deuGenfer, name = "martha")

# Polnisch, Indogermanisch
morph(pol, deuGenfer, name = "martha")

# Ewe, Niger-Kongo
morph(ewe, deuGenfer, name = "martha")

#Telecingo Nahuatl, uto-atztekisch
morph(nhg, deuGenfer, name = "martha")


# Gegenprobe Deutsch
morph(deu2)
morph(deulutherAlt)
morph(deulutherAlt, name = "david")
morph(deuluther)
morph(deuluther, name = "david")

# Höchste Werte
sort(results, decreasing = T)[1:10]






