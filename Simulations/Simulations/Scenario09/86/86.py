#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.85660173262268]),
array([33.221823886647975, 34.34989971203293, 34.92682969458409]),
array([32.10923663979911, 29.11144393361879, 28.977383046463515, 26.364074678640666]),
array([31.100311458644892, 32.51938230387627]),
array([31.95796193296074, 30.863130015529556, 31.4913354296642]),
array([30.93255027406359, 30.81453423775721, 28.91706330552353, 30.938906371068168, 29.81913454816026, 29.007157527808413, 28.01417086084282]),
array([29.664122137030944, 28.391066003791327, 26.514028573024753, 26.535857428313196, 26.22287559228671, 20.62634106866068, 19.37202997968273, 20.394983396395258, 18.792019315299953, 18.557227175903332, 18.94677946558038, 20.114152165380936]),
array([29.868331415917055]),
array([28.279677363905016, 17.30141723223034, 19.496671010384915, 17.237614899659615, 17.23195810465691, 19.7901614576703]),
array([26.877716502261123, 26.986683972072637]),
array([27.90322808538881, 26.946568495062113, 25.703519452999863]),
array([26.562981148909675, 20.748049013209446, 20.253280069204056, 20.301492268077535, 20.373939819227346]),
array([27.533725680623007]),
array([25.608021959168234, 26.372242812671228, 26.82838762397022]),
array([23.996778800892756, 27.23775880347581, 24.4768754850227]),
array([27.588946357648688]),
array([25.640275062438114]),
array([22.958834847948737]),
array([26.072395448111173, 26.30256552960806, 26.05923977343427]),
array([25.508784267491226, 25.59371212904521, 24.100809589125404, 25.274086279495712, 25.201357556720975, 24.634140970323706, 24.51545659391564, 22.93915092679314]),
array([25.66935330833918, 25.85120718367929]),
array([22.724370382266724, 22.8208809075592, 19.16392488046904, 16.352240476784733, 17.71263475616466, 20.13466892193709, 16.40726711713478, 14.520918066387647, 14.765001365032575]),
array([25.193183090505297, 23.84144245407376]),
array([24.693568443483173, 20.610230429582668, 21.667238709185977, 20.96340181319286, 20.370707647438902]),
array([23.90712530219013, 24.69999037375705, 22.898742026717514]),
array([22.93380602873046, 22.219214780323128, 19.814503870699856]),
array([22.84494317463316]),
array([23.683670066879497, 23.32130290186799, 22.524225345785776, 20.56118788968079, 22.69255072470119, 17.086804075516813, 16.261278165759393, 18.636540365621055, 18.18014283543568, 16.266860353120585, 19.58528477876701, 18.63658349880411, 19.93855551439436, 16.026211446905585, 19.35440331029882, 16.3087939390564, 16.96347001644824, 17.73038139702728, 19.59665712308187, 16.705236536380397, 17.400712446939927, 17.082431279979076, 18.03698182852283, 18.336672168338296, 18.80731621919664, 20.156962475033072, 16.778238920551193, 18.922089160044415, 20.111408562926, 17.976465490744506, 17.38425347633773, 19.94194446019412, 16.963424126983185, 18.337534035750256, 17.31134499788692, 18.22850994899649, 19.024513588800055, 19.42887936976329, 19.524724322649877, 18.41487961517901, 17.79360908731562, 19.936038765339234, 16.412381244727115, 18.721003168685385, 16.869922709200626, 19.793390803882176, 19.30707169055331, 17.293238856820263, 15.502277045969546, 15.189516647054981]),
array([21.118440617075983, 22.900482841967023, 17.81193129350144, 18.22002149527546, 18.41932124615125, 18.96338727858233, 17.79248711882162, 20.12475784289869, 17.50856506706929, 19.42828745811315, 16.849377153816796, 17.555267260355738, 20.40765210108036, 16.002447196616107, 16.377841668883114, 18.844127785971878, 18.699844593720076, 19.775324747856335, 17.09450056335048, 17.46262033684154, 17.145325521160586]),
array([20.485291092787712, 22.213818338036646, 20.244830645671655, 20.058955969638426]),
array([19.582419794577174, 19.683256192435255]),
array([21.8434409936185, 20.69133400718764]),
array([22.62484258403263]),
array([21.54813425822226, 19.981055525746562, 20.217437011917724, 16.566428749826574, 18.286370874334885, 17.935681091299184, 15.991882093130009, 18.652717115712086, 15.594906047418144, 14.7938188645613, 13.704509363548842, 13.051417172652437]),
array([20.67101351635902, 20.794676474537212, 21.631820195767407, 19.709904611567698, 19.76783478521458, 19.298026806754375]),
array([21.732196202493856, 18.96683185787673, 16.266870984232895, 16.45395348266558, 18.659312947399986, 7.264670278955253, 8.307391590485762, 9.946765628370375]),
array([20.489821782445347, 18.113594634315817, 17.070550094764805, 16.82584106370497, 18.653403633397488, 15.709942566688852, 15.051000856664377, 13.688224550141435]),
array([20.48423982210771, 19.66436318509913, 18.710450843117687, 17.488452713220763, 17.240598797560807, 16.396699090210184, 19.730547584881315, 19.527930381159216, 20.158184046638496, 17.17478556432039, 18.068920274790017, 16.41120747642035, 18.116336571821055, 19.48707538012992, 14.572488268630984, 14.076156588785445, 15.151537766570897, 15.38213472320804, 12.265933494811259, 12.827294047044838, 10.768583758836986, 10.94025001582074, 9.115063801487006, 11.557567147652124, 9.088164545800096, 9.897473241005176, 9.806019068196218, 9.159879189590933, 11.283217853608726]),
array([19.027248070509568, 17.114781006475756]),
array([20.473771523206135, 16.779957718511696, 19.948413310398802, 18.955281049781785, 19.331985403675855, 16.489058614996104, 16.10265885087784, 16.67452032305362]),
array([18.59885678050029, 18.522605050667565, 19.07879583240675, 17.656860551188526, 17.462131816740506, 16.61054144508107, 16.60170244618028, 19.281426794777865, 17.619772866937872, 17.300345575628448, 17.530525100659638, 17.296086497582582, 18.662492871081334, 13.992809579765444, 12.450013706423535, 11.645254286288463, 13.053846160937535, 11.728986535681411, 11.420007888114526, 10.668328979135978, 9.105714026805062, 7.7020528017717425, 11.136200717430878, 8.758431635979466, 7.274448398765827, 8.356326098723608, 11.00215970935739, 9.872997816653125, 10.362834683101504, 6.768269887098452, 5.858585341335589, 5.6573469116466235, 6.447608159931233, 4.728525630797132, 5.195058305986941, 4.161614190798119, 3.797895490534601, 4.032609928818344, 3.9568613177059255, 4.831947063406612, 4.8131246331475115, 4.413952975221779, 3.748582848766897, 3.617539373932853, 3.276257685502231, 2.617088724187748, 1.1566340283738985, 0.9610322098648254, 1.1052230452043847, 0.9771527783042746, 0.13217074113925165, 0.2955026828781466, 0.44308131635792614, 0.6422740368941999, 0.0]),
array([19.198738383310182, 18.66245308380783, 19.231933737963185, 20.1595292411677]),
array([18.50076859094269, 16.21272340848913, 19.852952903996446, 19.784032904547107, 17.19465149266564, 18.133789770682117, 17.906877699089847, 19.21924166296185, 16.165876278768543, 18.888877590040945, 17.646956857449286, 18.34413229906004, 18.53169542371937, 16.752361959840947]),
array([18.81626804158388, 15.998350392031224, 18.05778231360371, 16.30596227991502, 18.77845259314553, 18.97520215901454, 18.520384025555284, 19.213442391567074, 19.603071402669233, 18.42611106644477, 20.051312972052752, 18.326881379780172, 19.302674156117103, 16.096534595562197, 19.93292499234062, 14.771854545647603, 15.432870363427186, 15.054157246684792, 14.379624412925354, 13.73419741043598, 13.686701189723458, 13.428173580267474]),
array([19.76817322225684, 17.750946639035973, 17.693938784236614, 17.566501709548056, 19.90885510329433, 17.2085752946637, 19.438204942930188, 15.988435330707299, 17.807920088551718, 16.49691140776668, 18.105816951199728, 19.389612462033483, 18.188794589413085, 18.474219389300792, 17.504006583299912, 19.257399024227325, 17.21952959765163, 16.068533786499355, 18.368833936533747, 16.83820913992018, 15.770268964798646, 12.86369687057796, 13.755529402225715, 11.66043266469067, 13.00659307699129, 11.973078412859403, 8.109232508383265, 9.865906371101373, 10.470857575155552, 8.254706493971094, 10.624784594046847, 8.14652624252782, 9.123783850869248, 11.512593235051346, 9.99121206149675, 8.474246586355259, 8.022348096589653, 10.119962426650526]),
array([18.05492388166809, 17.966029042152822, 19.14096773260176, 18.994810961035004, 19.305383639635497, 19.472366736962403, 19.61964063523907, 19.108937740721068, 17.918395037509516, 17.948552022385197, 19.636582639453426, 18.613883627010654]),
array([16.10886697842071, 18.784043651856276, 17.681098361148724, 16.430911079538596, 16.57988134207286, 17.765074108953016, 17.462078897404968, 17.599611890070026, 18.961698465233898, 16.980288688998055, 19.11166188574042, 19.404237100170615, 17.153141513754115, 14.566694598175543, 14.199254384738602]),
array([19.028608331979097, 17.92726980979433, 17.76526872120065, 16.155673784986575, 16.535716713615763, 18.544121895183487, 17.47395867023966, 18.162476223850128]),
array([18.09354742570443, 18.32229759033287, 17.769371829053654]),
array([16.490057945725457, 16.176881529027796, 14.26598623534735, 13.434903270027837, 12.538001997915993, 9.43876003052134, 10.305495098738206, 10.92111802432429, 4.580684535588949, 3.5737597526857385]),
array([17.80037930162292, 16.21381615912839, 16.269130687016307, 16.625681024242883, 17.380968023541012]),
array([16.137665574526135, 16.174069276414507, 17.225247475990958, 16.772973514858283, 16.49232351128212, 17.405583212166896, 17.343184196542186, 16.35343895103592, 16.78306854895068, 17.66997621804233, 17.64468002329224, 18.11283197552398, 14.878331628251766, 14.28284900584851, 15.698365435029128]),
array([17.135237583815282, 17.609502139938694]),
array([13.316896096199061]),
array([16.872584872976283, 17.17361581811127, 16.999362409308905, 17.221879504201635]),
array([16.90884951317632, 16.832850647733817, 16.01966345009753, 14.40843163224953, 14.761565845728512, 14.52644855307674, 13.434711033936582, 11.295013189782273, 11.252736900987653, 11.431291567279423]),
array([16.750089928663577, 16.364719933713992, 16.02600209925471, 16.360104602586315, 15.555468016250362, 14.864789781843726, 14.31120476379572, 13.68820534552619]),
array([16.695091691525636, 16.18702938828046, 15.949006199819738, 13.778364838926551]),
array([15.971221628662356]),
array([15.53473048262463, 15.084918609480543, 14.552389835173734]),
array([16.03979221701745, 15.416678665306174, 15.390606882030173]),
array([12.886467365170294, 11.514031335782214]),
array([15.031935402505095, 15.285173209699428]),
array([15.101400874344044, 13.052990621532558, 12.003268686458517, 11.642514457096787, 13.173594956336565, 8.014676713330728, 7.621813878602568]),
array([15.190198980960732, 13.134433430090613, 13.20632230331548, 11.49419720524658, 10.39791155563052]),
array([13.679663727187918, 12.917558749460351, 10.212365880826145, 11.588201381005517]),
array([13.386405632294897, 12.499921726586123, 12.689443115134914, 12.160726559185452, 12.861581757844277, 12.795106031216996, 7.493590019821816, 8.657042189734899, 11.167459758483883, 6.535504545824761]),
array([14.60693037543189, 14.58910675296034, 13.746757372278275]),
array([14.271653316358762, 14.00588591434992, 13.779303585499932, 13.596341442220535, 12.417484901491985, 13.287808044697837, 12.356408694337953, 12.05130241000486, 11.91182561749026, 13.45812128388919, 13.244939035633218, 10.13552042709833, 9.958523932251385, 9.52868874256561, 10.270175676429885, 9.941785919218399, 10.714077975901175, 11.256852789627647, 8.982334382185282, 9.210081729683289, 9.844967192590442, 9.285796701105072, 9.626927770583958, 10.860233951815044, 11.131160202919935]),
array([9.879275831505009]),
array([14.299942574480271, 11.833029555228816, 12.08253284702011, 13.186829193120452, 12.527069784334518, 13.086211346982637, 12.386952150958782, 12.674930353022724, 9.84208701269695, 11.603376544165334, 9.78165995750906, 7.833811518709942, 8.109138805443902, 9.341844525372712, 8.218493655512475, 10.106771603373184, 10.798606970235587, 11.029363461834524, 7.065738852697162, 7.00707514875879]),
array([13.958377445021084, 12.094167409814284, 12.734585450574414, 13.603236527583174, 11.511754586657785, 10.999653597533833, 11.444414141358308]),
array([12.82246271652625, 13.107509159084488]),
array([13.938321033311475, 14.035873223711436]),
array([12.234507618672316, 13.371302907586404, 12.829149349828612, 13.30410001262358]),
array([13.150005616513713, 11.860778295411524, 12.721388101271824, 10.458667532883279, 10.682531873193227, 8.73303369229132, 7.446149278271632]),
array([13.682982121345672]),
array([13.894734979162017, 12.468314311116247, 12.284598556809382, 13.536848512763063, 12.048592989189517, 12.333628316026056, 13.64597536014649, 12.777789657691523, 12.168740889264063, 13.282778241573217, 10.212473300510371, 8.038516905658824, 10.282991800796138, 10.481813581825769, 8.449998517216569, 10.186499176540046, 10.605643689135476, 6.625176330622499]),
array([13.46783381725763, 12.309300171499359, 12.054931259866741, 12.524896328386044, 13.170732099102272, 11.690652568566966, 9.30561668785942, 9.784876434397649]),
array([12.756414675765098, 12.008880956665752, 13.470835492085762, 11.210560450406685, 10.548333622234042, 11.599733993617441]),
array([12.477945282066509, 13.19584845768909, 12.98160429151397, 12.581633267234322, 12.125777595439246]),
array([11.595940875213772, 10.318934675503451]),
array([12.653939808127642]),
array([13.041441441278158, 12.257743865183421, 10.904939293020522]),
array([11.92814926300369, 11.70266691254672, 12.876306503211879, 8.000475147709762, 8.377972024144356, 7.998047656009241, 7.436226119651497, 9.400218356083194, 9.39391374499074, 8.584822546915937, 4.784888520003966, 4.3374144090204005, 3.4526821548998]),
array([12.14531348399372, 12.209318478449731, 12.243029406897111, 12.26601714536031, 8.680893874984083, 11.170527172191315, 11.496953074392884, 11.053680504221898, 7.318267839627848, 9.509085325523083, 9.716240100281663, 10.862217139236511, 8.333804444341517, 8.445764030309327, 8.758963558248832, 7.210255593652334, 6.6898687563928085, 6.874548626628911, 6.098921698590827]),
array([12.051288792685806, 11.785388636086996, 11.324591255998566]),
array([12.496405755031109, 12.214459585831058, 12.378583327056413, 12.111247404530719, 10.909755761006938, 10.622204834336634, 9.709094301324736, 10.822264682890118, 10.386688371634175, 11.358703986287804, 10.380515280295866]),
array([12.070165244260838, 12.62185824204323, 12.455519949033564]),
array([12.022200983821923, 12.028604752798072, 11.63351681876466, 11.413546992853302]),
array([9.976989167589366, 8.442342564530945, 9.641522610147357]),
array([10.142530858004971, 11.534813082414257, 9.463710288453537, 9.774370083131501, 11.402304350637117, 11.620799261667582, 11.349307241402359]),
array([9.650487880473591, 11.618287481705217, 10.366524126444695, 11.308769238085741, 9.829780547384408, 10.328204727403698]),
array([11.653479515682601]),
array([11.002495175154335, 10.115166731082637, 10.384943188631265, 10.002671303025076]),
array([11.677349705399271, 10.093612069928994, 7.551753346888142, 8.365680333770847, 9.992802929242774, 7.982424704540577]),
array([9.994053084463438, 8.372518627366368, 8.4358569771999, 8.670728023284227, 11.451646065929872, 8.424703753818887, 6.09367776603362, 3.729563474747236, 3.8404727877618785, 4.776266069982126, 3.8966492991337836, 4.6991123132344885, 4.521937354326831, 4.58296520985465, 3.38348997349525, 2.762804972097068, 2.5666499429646525, 1.474580969418416, 1.3123607658275247, 1.4722148292757602, 1.6179356668353335, 1.2080127305464714, 1.7775467433543195, 0.4995325365799223, 0.0]),
array([10.821076342177898, 9.928595245087308, 9.534639185173544, 8.733140828292711, 8.004382865700377, 9.27944681074441, 8.489600969846434, 9.52706349137042, 8.413442620848459, 7.8057981381151915, 5.533796342863591, 7.148138418725824, 7.193475107477874, 6.2133029508448665, 4.698774130503934, 4.77200699265043, 4.766961761877624, 5.04775626980869, 4.979356589404637, 5.199276481950305, 4.751313518867969]),
array([9.31936381069019, 9.425186464383923, 9.473476952531351, 10.353509899817645, 5.853532187812771, 6.161179126638753, 5.30033816270916, 4.066400202441818, 4.903898749716711]),
array([9.70317346918142, 8.039356905783523, 8.529476550924507, 10.662203938746874, 7.422718631531177, 8.038260984366564, 7.134673189848193]),
array([10.412781409645987, 9.215911318483434, 7.359702100866114, 10.125279915817227, 8.304538145937293, 9.209922622941038, 7.9585336495338534, 10.00985616683067, 8.171820320739396, 10.623618950169368, 10.140978566591532, 9.973787867225377, 5.800130079974037, 5.724262119589492, 5.383180543120023, 5.728656121108251, 6.604329576314155, 6.0946255219393395, 5.689842204302199, 4.203050229436553, 5.317041726775581, 3.6664052773725224, 3.6329638434553178, 4.005567123917242, 3.9697800081330357, 4.525622729262975, 3.119186207137799, 2.909213934407403, 1.8367061775506421, 1.745988375914747, 1.5322008238461997, 1.7472789088931886, 0.9173429388393768, 1.4651333781789406, 0.41527581797383856, 0.7102864412361746, 0.7401776270598827, 0.19749779050666116, 0.7474187631470569, 0.7739890174717882, 0.0]),
array([9.825125131212879, 8.191484372257854, 10.54063726257446, 10.488083522726441]),
array([9.781250242355512, 9.055990673600007, 8.564817766198736, 6.967847751935854, 7.188004258568189, 5.675498665552577, 4.205430074978795, 4.152903069332211, 4.405116689645902, 5.026619413440815, 3.2833733127162197, 3.289679046285143, 1.4680053565620939, 0.0]),
array([9.634136660139848, 8.617554572907945, 9.988363586326539, 8.179478574888451, 8.800262650481779, 7.778881441333583, 9.769089495297344, 8.526665021971333]),
array([9.287766466363728, 8.410227376699718, 9.778391863120918, 10.09562721468012]),
array([8.054701293991707, 8.252304248169022, 8.007121685401732, 8.817341867673136, 8.30503912328861]),
array([8.14519325055554, 9.612285602925732, 9.341658481147551]),
array([9.09541579696584, 8.986465671776205, 9.094123419645214]),
array([9.578098485692882, 9.410707433516041]),
array([9.102567417960737, 8.21458906828075, 9.46405130055186, 6.1921331376972075, 5.096169278639639, 5.018467838635069, 5.10604645040356, 4.798823130371522, 5.300092882671898, 4.589377574710396, 1.6109997333353197, 1.4137650821510426, 1.4839095340808268, 0.7903758316571934, 0.5633142250197447, 0.09691626979812754, 0.0]),
array([8.592024355981213, 9.01470061611937, 8.481469172701017, 9.299090243552165, 7.827540962553117, 9.061457015937929, 7.4065686409920595, 7.354074308677696, 9.219346639613354, 7.558960824180565, 8.999934026142608, 6.7328646370271885, 6.802463513891958, 6.1931936068311195, 5.38483990291794, 4.896618607046831, 4.533272757269816, 4.629889898452269, 5.096502094297908, 3.707806611323942, 4.804599066959018, 3.016274327768866]),
array([9.172776550356845, 9.382965329453393, 9.163759320901528]),
array([8.951479670695333, 8.93205052537017, 8.28035008494893, 9.033811736541779]),
array([7.76996059680468, 8.06492608438077, 6.933166501938498]),
array([8.409312956901008, 8.364401377834943, 9.07258098459922, 7.933436759913596, 8.08705210672138, 8.8952478251926, 7.871976189146038, 7.998656667659511, 5.582622413426213, 6.560491606356964, 6.927506012871684, 6.8804324852197105, 5.854618580537914, 5.154502787040041, 3.7185639962555257, 4.282157268263327, 4.452044107498502, 4.8485027466572745, 4.690744323278705, 3.7702486799721173, 3.662106025285377]),
array([8.504921846915332, 8.002260774416463, 8.66780093720425, 6.853561733505581, 5.866062695038378, 6.7445371273258505, 5.867852755299787, 5.519527242635499, 4.969082445027177, 3.965294405596536, 5.099409583297563, 3.6035688628400053, 5.10023859290975, 3.652909230849055, 4.6535629074203895, 3.0961003887935656, 1.2268593431483392, 0.9791280870712027, 0.9897314079804915, 1.0545188629099798, 0.4246378214010386, 0.0]),
array([8.35414646485656, 8.728857877886556, 8.408216094548624, 8.970433621749747, 7.409514607622081, 8.212777850960547, 7.4196714499747145, 6.6411574608677375, 5.552978260545915, 6.32132350775001, 4.812853830552851, 4.5684299208896455, 4.608224024784223, 4.943230929559561, 4.892292770400462, 4.1707930894854774, 5.263138913767881, 3.932947221347494, 2.8809700629148467, 2.499961685243562, 1.9626961953403335, 1.0460082203712755, 1.0311048713964905, 1.7183030023489112, 1.3318316847811729, 0.3988271075110554, 0.0]),
array([7.4783065190409745, 7.273643436813709, 7.400419877511887, 7.221361234938694]),
array([8.053120319403988, 7.982912160359824, 8.126675954810404, 8.484226352104693, 5.447105747522439, 7.156184949467597, 6.550902643333012, 4.9035341277678, 4.6365265307106664, 4.427001195106692, 4.9817344062069076]),
array([8.380107465856252, 8.46272569293549, 7.58755938087447, 8.20748132549028, 7.637156492378944, 7.514361778296416, 5.826858695604952, 7.098924881784846, 6.8866562777367, 5.994538230083645, 3.624676579684696, 4.458376299802081, 4.888952888111492, 4.861452859615101, 3.6920295649113104, 4.590074205806415, 4.45343970399266, 4.06090959313259, 4.264795568087349, 4.868154893722495, 3.620721643419163, 4.939105112979528, 4.585892261201036, 3.807381601579192, 5.257435392074172, 4.918535552030955, 3.4713794576832213, 3.417892080030187, 3.47655250797848, 2.115566808994546, 1.56716793934425, 1.5519086708258951, 1.131128676927065, 0.8425296864813046, 1.7005593569607287, 0.6544275160681254, 0.7484100480485951, 0.3502592073326537, 0.512150529077809, 0.03691309111803885, 0.0]),
array([7.56244981052165, 7.737293801341574, 8.145636877470995, 6.24024396792256, 4.775539812853616, 5.152198896863385]),
array([7.293044137321393, 8.076283898117168, 7.761624116967084, 8.18106420537376, 7.912475640227882]),
array([7.760942730498737, 4.2606804861499565, 3.861514049711878, 4.485976444999494, 4.675088864954147, 1.558040327270888, 1.4536886108736873, 0.0]),
array([7.406786906749062, 5.279574246583774, 4.041014684700736, 3.4959217117501673, 2.628977702302404, 2.66375794455183, 0.0]),
array([7.648728719467459]),
array([7.972573168781857, 7.51878390208485]),
array([7.7708497358423765, 7.246103439024932, 7.5368029065026105, 5.714390698688792, 6.516860376222794, 4.857949302290516, 4.3889272635907135, 4.630616807301662, 4.657062727884763, 4.938980883689779, 5.1661007407730235, 3.469555243294441, 3.296379982229716, 3.5158458590907156, 1.5562724247835367, 1.3745350794870426, 1.0584530114502386, 0.1987926343423928, 0.3056234608329153, 0.4787080840823032, 0.6611481763931161, 0.6962960898333996, 0.0]),
array([5.896128938248555, 3.731253283814116, 0.8844724305823132, 1.0888399160854503, 0.0]),
array([7.653096968995721, 7.454898967639218, 7.6591375297968565, 7.2339183122296165, 6.9584512555604565, 6.191613731508081, 7.036940725633848, 6.916197658215119, 6.705049129567857, 6.35950444943929]),
array([7.494970746944058]),
array([5.514121790921621, 0.0]),
array([7.283409948044907, 6.3022980009088885, 6.41594360245956, 5.90314904519505]),
array([5.005686428560809, 3.9513155358829493, 2.1065109498146377, 1.0748220039309184, 0.2518301438728916, 0.5309041682714963, 0.0]),
array([6.289693564848104, 6.2572980557410585]),
array([4.410258420788859, 3.925496922960858, 4.194735693367509, 4.962719638930482, 4.742058510745178, 5.113655980775026, 3.6318379965696845, 3.9090265622581075, 4.014197826606783, 4.049723420735175, 3.3014152687571974]),
array([6.966716772155961, 7.0805814524542665]),
array([6.081335985074011]),
array([5.9645275583480375, 6.359421557441207, 6.790293232663095, 6.378743338246013]),
array([4.510891281906089, 4.87901462232941, 3.602159433257322, 5.10687256344648, 2.877027589313783, 1.6506172537489685]),
array([5.798935246511277, 6.229541329114974]),
array([4.250329103495381, 5.050107082541956]),
array([5.90904279309572, 4.255477362219785, 3.772053823461947, 4.88641283741364, 5.056264785806241, 3.220310852437067, 3.1989443720125004, 3.3005893538830575]),
array([4.1761657777601355, 5.232693069264296, 5.320318466618489, 5.316872837864622, 2.8686807334048883, 3.264262068222789, 1.3916038134459554, 1.2916284973507317, 0.9941033225289745, 0.7572482341147181, 0.0]),
array([4.76298583747989, 5.13020841419738, 4.669914200966273]),
array([5.119215426379413, 4.640649392899204, 5.015559154323125, 4.954421475985365, 5.097461673560251, 3.760443615830477]),
array([5.393770442501418]),
array([5.010919710638279, 5.141006345330423, 4.848757951776167]),
array([4.927283455236184, 4.880940922253517]),
array([5.194584302511488, 4.493359782906352, 5.01782974600965, 5.116097918180003, 4.811567579594226, 4.658702589757836, 4.787054192427092]),
array([3.9383074471983357, 5.107647595730126, 4.32337789016463, 5.227802822674298, 2.924766376029988, 2.7112215760937763, 3.577285188336238]),
array([4.330980812250232, 4.690818883656056, 4.150723005788822, 4.469553026930612, 0.9395381110516337, 1.3859307318388914, 0.2590263398824544, 0.45732100256854175, 0.11038805093870345, 0.04339086732923386, 0.0]),
array([4.65707344224111, 4.653204262556845, 4.678064538581334, 4.89429316663301, 5.089484823042488, 4.590612940062029, 4.665603984860282, 4.723240871924931, 4.635175720277233, 4.924462768561847, 5.092500079749565]),
array([3.6698429611277112, 4.012304206241781, 4.032193558799408, 3.9976655848446274, 4.406732735289499, 3.563897399219284, 3.475095447816181]),
array([2.811831273992074, 3.070842906587888]),
array([3.9148684386527357, 3.928371365787969, 3.8100999719833184, 2.30856343959867, 1.9993173457923623, 1.9697619540062805, 1.3114669027909787, 1.0879761997956803, 1.390495179092401, 1.571432257704635, 1.4003647217104762, 1.7645040826882346, 0.9803475308809515, 1.3237577518537944, 0.17364749491898512, 0.7724631964132236, 0.5037045845086126, 0.5851293342502675, 0.6823230583935589, 0.6412712138928628, 0.04881878542753185, 0.10970396462957455, 0.0]),
array([3.7598551300919625, 4.059679419567525, 3.6347253118183627, 3.732989635380968, 1.7574314573283292]),
array([3.7657687541270866, 2.7380623562567568, 2.5857809708686412, 3.240125591957399, 2.992276064188057, 3.427691359416, 3.136411955159873, 3.3893835804393695, 3.055235248808034, 2.370581500260996, 2.002430623814557, 2.0812831472331568, 1.1675285747559059, 0.9808279031097967, 1.0838234105483147, 1.7735015469275186, 1.490533708064097, 1.671943355629121, 0.43194070988536326, 0.5375826093284233, 0.6029167406508698, 0.2519992848905136, 0.26699724424663107, 0.0]),
array([3.315312097173601, 3.549996153007381, 3.5754935801571475, 2.910315946410378, 1.7413043266292862, 1.701367005029014]),
array([3.065217665294843]),
array([3.4468763989769817, 3.2538696945921117, 3.307522060920292, 2.9018140476049332]),
array([2.2457377491841304, 2.2646003392984837, 1.777743292293073, 0.803743034914816, 1.1637461812586278, 1.6629596538998548, 1.4019264285529949, 0.0432494117417139, 0.0]),
array([2.613833257801203, 3.314924909969517, 3.0226746289134065, 2.3552851261798464, 1.715775338925995, 1.212349123469469, 0.8128431802681207, 1.3211834712186623, 1.730407446385356, 0.0]),
array([1.4416013818677835, 1.2625508444611757, 0.49614556523336567, 0.24083444549798938, 0.0621907752909371, 0.0]),
array([2.6993898457181578, 1.486509723101829, 1.1023367956817989, 1.4963597891278608, 0.47345513201501105, 0.6519422448676373, 0.2217916010685853, 0.688905602376153, 0.021584469347956603, 0.0]),
array([1.554917480464454, 1.6550882299920506, 0.2946545209365952, 0.6148115411806683, 0.07294428143291362, 0.0]),
array([2.6417780837741955, 1.7113752222150262, 0.6531451018085055, 0.14208011681667376, 0.6443379543764856, 0.47180412534227606, 0.0]),
array([2.0287725288020813, 1.9385735108321396, 1.5138696607259874, 1.6234456692802044, 1.1584105889095497, 0.4357036197507242, 0.3269701293060193, 0.7012376555470425, 0.0471848883267742, 0.0]),
array([2.6057268838012724, 2.2018272224075526, 1.7942932556355797, 1.6664617848388426, 1.7179974812000909, 0.23944573346847864, 0.0]),
array([2.0080946005166487, 0.0]),
array([1.6031102844940448, 1.4903480921722752, 0.3076721209008933, 0.6077659656592895, 0.007909056049249644, 0.0]),
array([2.1973075365008987, 1.4643835256218203, 1.662072918496408, 1.7549146661453667, 0.5785317733407905, 0.6405151659315136, 0.0]),
array([2.148712504985367, 1.4135883701694596, 0.7095693470550153, 0.7341606166049031, 0.04600227679413531, 0.0]),
array([0.9075241939366975, 1.6290833993789051, 1.5164645469631544, 1.7542385664033804, 1.3328901941533675, 0.5812860755088851, 0.0]),
array([1.011995924529169, 1.703388769919557, 0.274793024016027, 0.33745985904909676, 0.39261859684325745, 0.7741813520937162, 0.030452785595680418, 0.0]),
array([1.1453395990085662, 1.6596431134458365, 0.8669837286137632, 1.7229926802591105, 1.4769096700102708, 1.3364468139025463, 1.018694266189423, 1.2177173136381536, 1.0571213445484415, 1.6515685698967315, 0.5662302897988939, 0.14484577330840487, 0.5737328757710394, 0.0]),
array([1.5885396518785027, 0.17297489868387916, 0.0]),
array([1.4373432824849506, 1.4836973478084161, 0.5320201205692692, 0.6914178617002417, 0.545829894386747, 0.0]),
array([1.2630533566204998, 0.22532946005160037, 0.09057639542902202, 0.0]),
array([0.9930510404713483]),
array([1.0296694044455768, 0.6945791836779492, 0.556368767419406, 0.08188642144048178, 0.043697168836797834, 0.0]),
array([0.42927447905967003, 0.7733042771665563, 0.725787172954049, 0.304799256611992, 0.29344229763478075, 0.421280473486578, 0.0]),
array([0.5781574893498423, 0.0]),
array([0.181657188543616, 0.2870622129633905, 0.05172812445494002, 0.018285327757512918, 0.0]),
array([0.23371509853589306, 0.003262797687343305, 0.0]),
array([0.2975535510760152, 0.0]),
array([0.20788288347956152, 0.0]),
array([0.004417518455955741, 0.0])
]
d = [data_1]
names = ["86"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T4', 'T9', 'T10', 'T12', 'T13', 'T15', 'T18', 'T19', 'T20', 'T21', 'T22', 'T24', 'T25', 'T26', 'T27', 'T29', 'T30', 'T32', 'T34', 'T35', 'T36', 'T37', 'T39', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T46', 'T48', 'T49', 'T51', 'T53', 'T54', 'T56', 'T59', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T68', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T79', 'T80', 'T81', 'T82', 'T83', 'T85', 'T86', 'T87', 'T88', 'T89', 'T90', 'T92', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T109', 'T112', 'T113', 'T114', 'T115', 'T116', 'T117', 'T119', 'T120', 'T121', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'T130', 'T131', 'T134', 'T135', 'T137', 'T138', 'T139', 'T140', 'T143', 'T144', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T158', 'T159', 'T160', 'T161', 'T163', 'T164', 'T165', 'T166', 'T167', 'T168', 'T169', 'T170', 'T171', 'T172', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T179', 'T181', 'T182', 'T183', 'T184', 'T185', 'T186', 'T187', 'T188', 'T190', 'T191', 'T192', 'T193', 'T194', 'T195', 'T196', 'T197', 'T198', 'T199', 'T200', 'T201', 'T202', 'T203', 'T204', 'T205', 'T206', 'T207', 'T208', 'T209', 'T210', 'T211', 'T212', 'T213', 'T214', 'T217', 'T218', 'T219', 'T220', 'T221', 'T224', 'T226', 'T228', 'T229', 'T231', 'T232', 'T233', 'T234', 'T235', 'T236', 'T237']
def get_taxa_names(): return taxa_names