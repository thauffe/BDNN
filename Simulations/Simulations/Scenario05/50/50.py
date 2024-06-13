#!/usr/bin/env python
from numpy import *
data_1 = [
array([22.577773530218387, 19.63413950569713, 13.691712896857174, 13.50582389785562, 12.128657948172902, 13.230391529739542, 13.042018591632152]),
array([31.762244518049005, 30.953738513912036, 29.3002522282997, 33.48027929248758, 31.345383858936373, 32.62643556410327, 31.695629542838482, 28.665120941764552, 28.292771707949942, 31.25268003850972, 29.207311333367638, 29.393856950723062, 32.65800901080878, 33.553051331169364, 30.17270205784204, 26.29194153494535, 28.072608523523314, 23.995312882062237, 27.222217474981136, 24.772500325101046, 24.141122660491384, 22.974832483470614, 22.440677015530692, 16.037463611822005, 17.49919161044617, 17.50175486881398, 18.335010146543848, 18.050789102423504, 18.04219598472955, 16.56164825118767, 18.24139562503059, 17.389593744714702, 19.910120279707325, 17.45448502547436, 15.70430002875969, 13.411038865226006, 12.737885855923746, 11.654320968101866, 13.394409898409105, 12.504155548237524, 12.40927004815401, 13.382331004312991, 13.370620207253813, 11.655689091443396, 12.262060237872573, 12.491621684624135, 11.767633164958642, 12.114039898909402, 10.97823015163013, 10.766669564572506, 7.764634011778783, 9.252772995292085, 8.200258541883976, 6.618914624809215, 5.564660624478181, 7.199822655412995, 3.88543059780927, 4.550260983786968, 5.036559267506946, 4.109375919655033, 4.687467891592124, 4.337687050054342, 4.526662766690273, 3.7787171658346193, 4.496786384320836, 4.827013395361658, 2.060390082174928, 2.465010424808499, 0.7983953379936359, 1.6237593594725839]),
array([31.520915819457397, 33.43296419990136, 29.77076970261432, 30.364172069955508, 27.0166649609276, 26.66275890827701, 23.509627707431974, 21.369418437993424, 22.13800724253775, 20.46922596933459, 21.947927841824097, 20.935248823906996, 16.40957325680582, 16.28763548571616, 16.535606185037622, 16.38778972888919, 17.832426782017716, 16.77001092834549, 15.730735605206972]),
array([31.32414085417859, 28.92966673353008, 31.767331489860133, 30.27160128323092, 30.535686043370415, 30.843815177611727, 30.459498534660263, 28.11177236755648, 28.503573749385204, 31.185455163746106, 28.296429073893535, 30.875121338484142, 32.114181497461615, 25.27777209853755, 27.087230355747828, 25.872740021325754, 25.21074807054483, 27.875772157759304, 27.767435677550324, 27.23674618629775, 25.59548113566498, 25.398352348039474]),
array([29.601319550305803, 30.54502517932937, 28.139211626541375, 30.300422087855935, 29.8148733400859, 24.77255146584229, 28.08964588761457, 22.95456548679644, 21.979308607716415, 17.032738055195036, 16.734876889243257, 17.41892393757893, 18.810435696398073, 15.121041806916951, 15.105587191998149, 15.403938346731895, 15.887358276028598, 14.685161445298865, 15.674515840092761, 13.086802534639268, 13.69248164736015, 12.84828117750591, 12.74038389239199, 13.548902643313674, 12.88831247729064, 11.645547063268772, 13.278304121757143, 12.007572915580822, 11.822601650452226, 11.830361948900135, 12.174627648732626, 11.756838082459232, 12.301157508785405, 13.714976933465532, 12.04192654097998, 10.938343539398875, 10.504707708600602, 7.808753611872914, 11.615402391594477, 9.931284956988481, 9.904546009187303, 7.897933802536258, 11.16973577170353, 5.785685914288807, 5.503623123532732, 6.414465929930586, 6.957144997058772, 6.94773015916786, 4.310016207457682, 4.382531725571846, 5.010341982391373, 4.552520665578307, 4.131394652502109, 3.7895805073476394, 4.394064516857119, 4.466147865132219, 4.58550722203972, 3.474162722967413, 2.435154809999686, 0.7871697651046328, 1.244361590778834, 0.24164769904397831, 0.41348579177714473, 0.04652547056716476, 0.0]),
array([29.168526519801777, 29.06334043263366, 29.739590088398394, 28.48069473808541, 29.557156454120936, 29.418137118546376, 28.389591879882772, 29.879698478326773, 29.750468756914902, 27.79643689178761, 23.83797261428019, 25.316371684102013, 26.323203283594705, 24.58817882579315, 24.636417642392885, 26.375467239404692, 21.08288867172169, 22.690488466935545, 21.939650468814005, 22.21240572422143, 22.45603125307738, 20.74285297505677, 22.54026014091649, 17.724117339340467, 19.33072107339254, 20.27800132010812, 20.36515116060909, 19.162271172993957, 19.821888959633824, 14.264830969270411, 14.25614904185219, 14.381165664004179, 15.55072901841843]),
array([24.79712118893025, 25.579138422590766, 23.371964446220638, 24.1727049286557, 27.707858734570188, 28.00161201867933, 20.440881109956614, 22.946180755258947, 21.258295072578363, 22.72617762580181, 19.837143960619713, 20.277672331553994, 19.840832554122134, 19.164545067936423, 16.769200788135358, 20.081397933318, 13.944682449681391, 14.853363527091535, 15.566712607714432, 15.812870840297325, 11.802273581101119, 12.586777677297714, 13.277748671379102, 13.470226330637194, 11.952511809972147, 13.469374568430084, 11.715649978918465, 11.806242770033142, 12.912062128697217, 12.455744816890606, 13.75061627475124, 11.799770309007721, 12.206339753855419, 13.02287628023562, 13.093149154603813, 13.135740427671633, 13.313145571068857, 11.4047320827464, 11.609379180857177]),
array([24.06485443364016, 21.610088845689035, 21.522838363201306, 22.393759100111076, 22.6049753419007, 21.421879246747196, 20.65314386755354, 22.17330602076135, 20.261282772294898, 19.11251206172613, 17.351467181962178, 18.509847252338062, 16.404870081597842, 19.18637255668837, 17.490503064871547, 18.196771665001236, 17.072231024597087, 17.729375645598314, 17.87064117979406, 17.49685145806229, 16.849261508933104, 16.96691054187624, 15.807120017837331, 15.468529341804992, 15.270380876273054, 14.558556096188921, 14.762282244682222, 13.841845892608003, 11.801927273222317, 12.440681551302422, 11.878399268514944, 12.318188639985982, 12.095768140336155, 13.16098451652403, 12.418052669137381, 13.534796699631668, 11.826237648591196, 12.392986741249116, 12.524439540062211, 12.360332292162923, 11.828369209708347, 13.614030189225446, 12.11865399767733, 13.225077013865832, 12.750564076773005, 12.758793674392722, 12.01839865213576, 13.157032982603983, 12.265487325061446, 12.153796320101552, 11.842265830867778, 10.764895704989751, 7.657366239779128, 8.410817667411028, 9.785753065181309, 9.97330266389385, 8.349447882505721, 10.67751451476428, 8.795903639036236, 10.720867386903846, 9.779384820256814, 7.446646544521622, 7.481508940455387, 10.336576651753706, 11.368723631120778, 6.9920028578278135, 6.201404972407133, 5.968520296446346, 6.628039585850597, 5.539599185363808, 5.781340032568564, 5.346343510949195, 5.735224831753194, 7.203036076815211, 5.338988302400578, 5.799053249835246, 5.627600478763597, 6.42631290380084, 4.968871261006361, 3.672210904203995, 4.8634559973534595, 4.944489176220445, 5.3159588349940226, 4.043173890024583, 5.210449918837537, 4.607139308209643, 4.923233396126757, 4.339994342159534, 4.615757136264751, 4.94139027172573, 3.62557775884774, 4.303917344808536, 5.210572170805562, 4.9409962262751845, 2.7704227535793873, 3.1504695453000116, 2.8280979463712552, 3.3471012577855124, 1.278140310975433, 1.2639336066138425, 1.3176425004820704, 0.11218747798992391, 0.044482151423245145, 0.05498348757422376, 0.0]),
array([24.49322471995315, 22.80221181498696, 21.14140981678884, 22.55905939656516, 19.597823245598825, 18.757963075119626, 17.937246142803254, 15.81222347065879, 15.556054714956836, 13.742739225458172, 11.757364575649586, 11.66065021210473, 12.208337007210334, 9.491674524531902, 9.210942478168946, 10.390686289308272, 11.25853393236181]),
array([23.93738876068762, 23.06469706121932, 23.97557133085015, 23.942426557653068, 22.281915641614393, 21.327733087047545, 21.49647978744385, 21.792903257063216]),
array([21.41297426029831, 20.392170962337925, 18.527494011971005]),
array([22.81523110569886, 19.059699962779472, 14.607087387758638, 13.494866353166543, 13.571403395666037, 13.30916017560314, 13.594792991346782, 12.982120290492935]),
array([21.1416003374791, 21.829803324661114]),
array([22.075412031243626, 21.739000313217275, 16.757813350738367, 20.20770501173974, 16.742472515621426, 20.250315861394494, 20.16420945218158, 20.415431163694645, 17.261612439526502, 20.37721347715024, 18.212385851637404]),
array([20.231110796032546, 20.372138531971352, 19.87099753051217, 14.439718205041908, 12.959408769787066, 12.828992349006962, 12.011090297876471, 12.064922967393017, 11.734655131572957, 11.682744107558994, 13.326007334590516, 12.753036259081444, 11.888047272467764, 8.523999610705138, 10.13889088036698, 7.702393811035607, 7.557294242777323, 9.69364740436296, 7.213348556718162, 6.891220515309762, 7.113072688543692]),
array([17.975396041062314, 16.225577981809334, 16.709835898323917, 20.079559598648423, 17.26573362250101, 14.190392832014869, 14.483014465752792]),
array([21.00915960686433, 20.962155407806396, 19.572122438838992, 18.800628713721117, 19.549328309351946, 19.107111788257736, 19.04523147032664, 19.055230734018565, 19.839996141635122, 19.9229152482938, 16.354653712198278, 17.331763390366724, 14.268215019145172, 14.904748635303747, 15.778036171721556, 15.106467371368627, 14.357850069698124, 11.82963506702077, 11.659031132388005, 13.472717214289167, 13.54971579414286, 12.483639732223462, 12.684726059173308, 11.892657590468833, 11.916865451319111, 11.828195868150814, 12.49545577345232, 12.028896269847847, 12.077442684715999, 12.180542081296421, 13.814935545861761, 13.337489747044502, 13.206761447465, 11.7065086456236, 12.75647625070896, 12.420785511138082, 13.754878864386598, 11.848696440708883, 12.962627760254204, 11.916725473263996, 11.86459200919854, 13.792617969904702, 7.612963970312092, 11.222699371060221, 9.095900299676838, 9.013324939655286, 10.706794174120468, 8.107100329542265, 7.765431046946318, 9.118955216200154, 10.251013709238919, 9.592907146957828, 8.872471756706526, 7.81616644212364, 7.929340069348973, 6.535709642480766, 5.348895714708246, 6.216738135040479, 5.809493435083387, 6.657811857249194, 6.1250971735771085, 6.404621622486153, 6.1000003377835625, 5.353447531397318, 6.564585429109407, 6.62016815653817, 6.776235184387074, 4.189067633496052, 3.7467536371153525, 4.28043318510905, 3.760786220925535, 4.895489584143725, 4.3990632776704075, 4.569144170811134, 5.219024202513477, 5.170887137381349, 5.0408453788200385]),
array([20.70186509612654, 16.6971976304111, 16.772157714493037, 20.371920478191203, 17.718500381059805, 19.5618911842213, 16.930737007902767, 18.007429641239476, 14.74451865309565, 14.562301735403748, 13.628986279198777, 13.08357911900995, 13.574820858517883, 13.53787135852535, 12.050604302977838, 13.015823172681419, 13.797831903896004, 7.58626820739503, 7.204031038702716, 6.313399643458666, 4.149273915623362, 4.044975568755241, 4.125151100550782, 4.299945984738763, 4.199114015866193, 5.102381064944876, 4.325536701881324, 4.615698206337497, 3.0281319398812467, 2.448974573935215, 2.023196156551429, 1.7066971649527285, 0.0]),
array([18.899987705739978, 19.29892786524197, 16.558743002337216, 16.550099246819755]),
array([20.72715081610868]),
array([16.595850231329322, 16.389191746580018, 18.102694879520552, 17.801446606228648, 16.099339252005336, 18.674417679280086, 16.50432369574766, 17.52501794079495, 19.197638721297388, 14.883646229104723, 14.898857613864386, 14.247439261110962, 13.94790629401422, 14.279484156853677, 15.495454892263565, 11.889570881785303, 12.516420507307426, 12.121128675869954, 13.607914829770504, 13.55604996797309, 12.150360718492456, 12.969827734224996, 13.092877351748164, 11.729328327208941, 11.65181515118079, 12.64620213816962, 13.05120411400977, 11.849329997798487, 12.413863494478797, 12.88644971169342, 11.688016511835649, 12.913206318370614, 12.208345108134814, 12.914848194078882, 12.78963606508863, 12.265411648806479, 11.646360488981234, 13.452604407619388, 13.101085022201088, 13.36204495489325, 13.720809457559039, 13.133484692454497, 12.600549008661202, 12.565772836878553, 11.763778722317836, 13.7022171953915, 12.410136545828049, 12.801528491261623, 12.27180545711649, 12.258150131450872, 11.605154844543456, 7.320959507495966, 7.6365178158076485, 7.563019926561701, 8.365305633402382, 10.13135656461938, 10.527021696827731, 8.864223905190261, 9.014954810961164, 8.778568003477877, 11.156210131399373, 10.775520362678678, 9.773134231569612, 9.552430278888368, 9.775021103685724, 5.973410563021454, 6.328017255854448, 6.722424157590914, 6.819664357908631, 7.209017517255429, 6.5124414154958945, 6.780980378848638, 5.619876440555075, 7.205213636531109, 5.499695204204743, 6.490445211411359, 6.605012927021927, 6.318218677215793, 6.008005425390582, 5.706394861290288, 6.1082948027774675, 6.631307559522739]),
array([18.283084351441882, 19.264877368303065, 16.531261462169415, 17.122844936266105, 18.779232526059534, 17.393693394218413, 18.724718904019724, 17.030314821807913]),
array([16.448934307024995, 18.63862473243302, 14.99903173190617, 12.71297640267386, 13.158046950622932, 13.621459115052248, 13.695754062118983, 13.563548877362798, 11.658757790224838, 8.354444289304153, 10.69070638896007, 7.957338656653076, 6.621811965470119, 5.94966525843363, 5.78126005348466, 6.5563622963628, 6.417145219801089, 6.270965342596424, 5.984503964979771]),
array([16.137024618211044, 16.81999581380267, 15.043856998677947, 15.348696593399199, 13.639515759807882, 13.080163892593452, 13.010116888197869, 13.439806428826495, 13.239188183450068, 12.952519034254896, 12.947760118676323]),
array([17.29401628799449, 16.796913196323874, 16.2568173324269, 15.197668137213466, 14.353817430069405, 14.29767321784363, 12.539900402606118, 13.466038467501852, 13.452174203947031, 12.007990784647365, 13.409894969934872, 12.494496610122189, 11.891589818576263, 12.985336929461452, 13.24553468599172, 12.477352975021677, 11.698108637426094, 12.266373472145883, 12.477337357133752, 13.146178825898957, 12.558270035564068, 12.50239854712844, 13.274009789212645, 12.821323538276857, 5.895343688411363, 6.928773146766799, 6.785665884195819, 5.800026146808401, 7.036675175877155, 5.557891555888231, 7.1686074872590995, 6.748738309968313, 6.935762785950421, 6.285556951773563, 4.499155852583846, 4.034428535757458, 4.847401757084463, 3.6974963377867662, 5.05224559892659, 5.26266597286236]),
array([17.00140585769774, 14.050928379992994, 11.767689669499475, 11.708872651965233, 11.950670139394735, 13.062914783248988, 10.378969411334152]),
array([16.029718189801226, 16.192129901191986, 16.851687925168054, 16.77642365997593, 17.387312592782894, 16.580605716813132, 16.85347969223643, 16.351605691939778, 15.35272425953537, 15.909836147621174, 15.885366843612298, 15.773055669118481, 13.914355374350743, 15.922781496870744, 13.125721137775418, 11.740810320957987, 12.73450275947129, 13.226335874062617, 12.717368695738331, 12.797210293432348, 13.807154219501744, 13.567761910486167, 11.750650174594892, 12.905598896635793, 13.617957981882531, 11.816077416398144, 12.126561361033971, 12.358224916397473, 9.805045339222355, 9.23110282825262, 10.475778949102873, 8.298287096279326, 8.282982411782386, 8.889807200439206, 9.155901804916482, 6.155029098266066, 6.700991417924262, 6.201987449976814, 7.094128118901063, 6.231834650426955, 5.6044787527415885, 6.598827703571793, 7.0292887829475, 5.884646269156631, 5.7565869297573915, 5.69467693465073, 5.025952936781015, 3.8562489920612695, 4.462586009371856, 5.324672831239829, 5.177491253820084, 4.189909009262436, 2.7079341771026813, 2.766698438602516, 1.9338703705342621, 1.8387930584005845, 1.4395165992564198, 0.864503861079976, 0.5883731399045826, 0.6464892619368844, 0.7727147259921302, 0.0]),
array([17.022690615700075, 16.141929644438235, 15.313071726850636, 14.259976469371274, 15.793476606701217, 14.368202523240194, 13.71274663001323, 13.619368093202135, 13.56778620031808]),
array([16.025135499865254, 16.096404415178345, 15.308899494688541, 12.465125818621383, 12.085933372039086, 12.745943784408107, 12.571462281251799, 11.851705335353724, 12.630635707469017, 13.601721676179173, 12.418159690643312, 13.45303442711045, 13.799898638901578, 13.090945071121245, 7.737958838766623, 8.339721612730708, 7.278487149435997, 8.919216027829847, 10.455944817673624, 11.142276159616475, 7.674358757755523, 7.735318591657986, 6.317048925362299, 6.233201219654734, 5.897544891655279, 5.791926832478361, 5.964847719289617, 4.895747615804565, 4.453147398498023, 3.633257908914892, 4.95433590823166, 3.3040096243602957, 2.9442557082317875, 3.0200742427103986, 2.531491833561209, 1.9500986614405382, 1.0912363169786827, 0.008578270345454816, 0.11353636806801631, 0.0]),
array([16.402297389231677, 16.203661218844022, 15.574015047299143]),
array([15.111172721710256]),
array([16.565714531832956]),
array([16.36834209079356, 14.447502681716847, 13.391328752599122, 13.175710845988242, 13.333876100590027, 12.405428973077095, 13.29933054960154, 12.774490268058052, 12.62516832801956, 12.45487227930206, 13.465106471572366, 12.306421866107831, 13.791349862488378, 12.77774281737243, 12.026354281520875, 12.896565537989023, 12.86619466511365, 12.25464932244795, 12.829643870798474, 13.31320552218852, 12.534340921901771]),
array([14.919626003382648, 14.40940954851071, 12.605499082025878, 11.805506971128954, 13.554468631845836, 12.466611701183004, 11.653322473987378, 12.403612353336955, 11.738481942458593, 12.914421417270686, 13.077184110205872, 13.09210652019625, 13.481976629956444, 9.661443096515022, 10.379615305698575]),
array([12.044822443433633, 12.794903195024572, 12.012350956015364, 12.189894193584305, 6.2549754141890315, 6.717000456205703, 2.8983252295641986, 0.0]),
array([14.146784657795028]),
array([13.98054118024689, 15.486784544405474, 14.771623284276076, 12.301437060166682, 13.66168640424936, 13.121399974529718, 11.658091361222578, 12.551926079206618, 11.677813338387004, 12.246579771859553, 13.572838361025022, 11.804050847919878, 12.828539210123093, 11.343770498052216, 10.355843832179634]),
array([12.103125606331197]),
array([14.534952336096628, 15.216173924161621, 14.23056874055342, 14.58994567786428, 13.860648037051424, 13.609176293099294]),
array([14.632268445988675]),
array([13.923218004322292, 15.183744812154048, 12.140427486706336, 11.856172138347873, 12.99212432370204, 11.664880673565644, 12.598961416815406, 12.162045658679853, 10.916962610437281]),
array([14.94898353607261, 14.087479398206202, 13.721119430191585, 12.154475556887467, 12.757144452410667, 11.73197146495004, 11.743635359855208, 12.759535358613046, 12.421558950047384, 12.003396065321011, 11.889382093104212, 12.5358678811421, 11.78930715373295, 13.743847399989805, 11.632215886700823, 12.468579043225535, 12.187428585364056, 12.248251984849912, 13.248283584267124, 10.778703318122394, 10.293654069322358]),
array([13.9867622321049, 15.39477214174843, 14.995836151886524, 14.463479559799504, 15.115700623012712, 14.592066315360139, 15.271432023846165, 13.834609952830927, 14.435404498255805, 14.469871441408808, 13.518734553854085, 13.104217669613686, 13.312038388270068]),
array([14.697079144424336, 12.976691856471323, 11.726854480351124, 13.758965452116675, 13.403235637892028, 13.247350744975789, 13.580011651072653, 13.57390608823879, 12.75463637087142, 13.17567810048432, 13.0670282971607, 13.515900903279947, 12.545064851074642, 12.529070248055806, 13.459505436367854, 13.225661414708869, 11.828961432348382, 8.527442526575495, 8.664108477195171, 8.746995995966053, 8.956240047012873, 9.784058227550029, 6.530749519977646, 5.361647667662208, 5.99726514465384, 6.605675950613559, 6.3371788472985795, 5.855971958526575, 5.621073639417107, 4.2152257933269865, 3.6663332545532636, 1.3929525484022691, 0.772331143229221, 0.0]),
array([14.628218550881401, 12.045735404572142, 12.494170775348376, 12.406184757389997, 12.996665453849273, 9.913169265110216, 10.678026948415868, 5.856686681134195, 6.56461664499978, 4.637202554277696, 4.432010043156953, 1.9496830302569537, 0.0]),
array([13.270750871577862]),
array([14.667791400474796, 14.279755101182877, 13.415302050766142, 12.762139128879461, 12.326069216134501, 12.252348459865612, 12.35899316112517, 13.102956044355421, 13.404861472461604, 11.85780291204856, 13.732651615614621, 11.73871104538259, 12.564287422853536, 12.516572639857477, 12.445841641720655, 12.19667391612148, 13.791220378253472, 13.45663517530188, 12.399862131121232, 12.602624487429733, 12.270982906469447, 11.950303562239645, 12.936164427982696, 13.508214301935343, 13.365342680265863, 12.378201994337457, 13.457552318098216, 13.400753235219618, 12.478724971202427, 12.513617414233169, 7.538617523400041, 10.366763662138283, 8.696496104834837, 10.44116169339893, 10.165702665795097, 11.251361543777044, 9.34446624898335, 8.351796126589111, 8.039481047125253, 8.550813631018105, 8.03755759068647, 7.491521357123695, 7.762149092884977, 10.63150050359498, 9.637415207322013, 6.613629476024404, 5.384120887840528, 6.661410349449219, 7.171672933310994, 5.96321268717023, 5.890339995198581, 6.993542943035388, 5.53425811410512, 6.453749404991016, 5.500504731016328, 6.301479017066783, 7.1540057110310125, 6.313251654335596, 5.286412912237769, 4.733046035418727, 4.399946401543752, 4.91031001682734, 4.3117506667576615, 4.735911089569997, 5.222184000309399, 4.706069467000363, 3.080079733584805, 3.1293166900812137, 3.19515681684406, 3.20525459031227, 2.059082307592905, 1.7151900490361427, 1.1192651131185762, 0.10226664891598901, 0.0]),
array([13.349196840510823, 13.205718966112709, 13.49594223048961]),
array([14.34160795976234, 13.38371700162035, 12.066601501906154, 13.467291321539394, 12.658876432849596, 13.290057952033164]),
array([13.84134899886405, 13.328742172075707]),
array([12.678054956065429, 12.802289642479797, 12.035941566152113, 12.715018032022732, 12.578880016541861, 13.091194525835913, 11.809373297659551, 8.835648929722376, 10.295687311599007, 7.403959293603459, 11.336200640308446, 8.615795865630638, 8.644359922150661, 7.525753234163254, 6.999617934556458, 6.385636050923362, 7.180619304712962, 5.757775127912154, 5.990152188524573, 7.020799935746586, 5.786488528902881, 5.741938826540125, 7.022101485628576, 5.6723540997298425, 6.685998179673058, 6.878859886513276, 5.256573534326345, 5.304897004511869, 5.309040133311611]),
array([12.963453893931908, 12.601666776457527, 13.48238151227909]),
array([12.610882066766797, 12.691248486994812, 12.137149683306307, 11.852794125886323, 13.05467745495099, 11.418335032273182, 11.624171020909978]),
array([11.967698475417425, 11.899625183978555, 12.228609125546432, 12.782011881711579, 12.259983761720711, 12.823787822675035, 13.092292396264282, 12.316274932646396, 12.109667477287392, 7.360907550013214, 10.272526382455704, 7.446581880135318, 11.43255423319608, 8.142115360798892, 8.209132015758476, 8.862211850666139, 10.308620281614735, 8.28127577980777, 7.498143626222717, 6.950131415817628, 7.029715208592576, 7.066401103593546, 6.855891809125482, 5.9953408832728226, 6.449126552549735, 7.157799826881293, 5.665067639241395, 5.780875347216277, 5.470325238204208, 4.464975319456675, 3.731774840696052, 3.935245105578275, 3.942799795220224, 4.217836653529174, 4.3037458759458085, 4.496298404606786, 4.7177559005332075, 4.180479480086095, 5.093028382990834, 3.8916854732687467, 5.090096685698998, 2.7482636703636283, 2.049853313328482, 2.1542759357420898, 0.823953335606574, 0.2230426069884064, 0.0]),
array([12.953766661348816]),
array([12.415394121978927, 12.259468653696231]),
array([11.803682411680184, 11.924795573679893, 10.16027790053764, 10.36091026373143, 11.0976416161652]),
array([11.708293109174727, 10.759222989479454]),
array([8.51369945164132, 7.834821449709686, 7.805712805089361]),
array([8.299330733850926, 11.12867568386644, 3.798243436339589]),
array([11.806557733613097, 11.831788605382568, 11.396945251557787]),
array([10.530059533007893, 10.93869847983479, 10.42536667901155, 11.06381048223152]),
array([10.902864476675528, 11.119006198682293]),
array([11.334603009813343, 11.426118147774552]),
array([9.558716487421808, 7.630860900334064, 11.11473948074001, 8.335158488482332, 10.06679883379539, 8.942250582741632, 10.900334746967616, 8.191649542103187, 6.784890072183814, 7.062071901704925, 6.878295577683189]),
array([7.788713095773922, 10.09286110239952, 10.320926597197078, 10.814646245268392, 9.392263004129449, 9.516655726338206, 8.25862382954247, 5.591233826373089, 6.634007338183023, 6.316951227411436, 5.388907799819692, 6.236997135535242, 6.384733869554219, 6.713930151080146, 5.8776162318712295, 6.463281268158083, 4.043153302909727, 4.844686428748523, 3.992576277195595, 4.854341300475969, 4.1162435812396865, 1.202142795740627, 1.7282367067331934, 0.30571807974056825, 0.6770065063566906, 0.00969430595043086, 0.0]),
array([9.277714596869526, 9.785648625197716, 8.743446787766196, 9.73562947403984, 10.154434702501375, 7.239287248749085]),
array([8.821030760938779, 10.703048430263252, 6.426247415423511, 7.048332157840422, 7.072615065395749, 5.952139143557112, 5.951005193121086, 7.054136165017433, 6.872175882456795, 5.1853932894176875, 5.194634547829787, 4.0968848782337]),
array([8.075460543705196, 6.496207362013186, 6.811064283258815, 4.521369723419454, 5.052215838300625, 3.9257840173170306, 4.958537569166503, 4.294487130834936, 3.798731332389622]),
array([9.190797967013756, 9.688440783605493, 9.891731106330948, 9.265387790861535]),
array([9.037493657079578, 10.692203580425112]),
array([7.487941423954482, 8.754702180574629, 7.9128209764459125, 7.027998225416235, 5.7281445096577785, 5.54255141809681, 6.872697452045743, 6.380710888412148, 6.045985080574393, 5.383629704572973, 5.455093531930603]),
array([9.465848721129772, 8.348387076788846, 6.44196929645388, 7.001387445922137, 4.719611216272714, 3.641627009924732, 5.208408298995059, 4.282103688245079, 1.998396305271529, 1.7858262158337315, 1.4198333553886948]),
array([9.443710816979104, 9.600823359418584, 9.265505968563593, 8.81144098347357, 9.186730024901912]),
array([6.914460353553344, 6.622949787533122]),
array([7.574381325417507, 9.109092772084782, 9.227908274588959, 6.491306742982367, 7.119127653564119]),
array([9.172881649967232, 7.33137862968751, 7.1228660268583495]),
array([7.375233914906631, 7.219336090597307, 7.2424131797662925]),
array([6.605031960608583, 6.727481798799803, 5.619567154272166, 6.996381298847523, 5.950093783172807, 5.100887883835537, 4.9943768544545595, 4.531232134981931, 4.53102899644709, 3.818734156433174, 3.773972658263463, 3.7383686147327886, 1.441379398919253, 0.011926694123812129, 0.0]),
array([8.694586162504509]),
array([8.452417745370646, 8.191903293666904, 7.086557743454568, 5.890954525602414, 6.115615244562531, 6.050901242755469, 6.660754980121056, 3.9326286511154445, 4.358391337718667, 4.8287519121742655, 5.295333130391538, 4.523403894194129, 5.3202989435726264, 4.989988608053206]),
array([7.419958697342361, 7.003216658159024]),
array([7.955543804269602, 5.9051007570780785, 5.8535600141703314, 6.464687392030687, 7.101799586035445]),
array([8.020964963897963, 7.498178632625571, 8.131667112566726]),
array([7.226978887592889]),
array([7.632597343028944]),
array([6.592729376617578, 6.020671807085512, 6.116013187855335, 5.627839412386565, 5.846105511814335, 6.832856513917105, 5.394737170059177, 4.214092740835182, 5.002148429751879, 4.210881750956659, 4.160845620792541, 4.628644014602421, 4.158322810915637, 3.9529019557469125, 3.754056119164853, 4.200855658174046]),
array([7.247638364635766]),
array([7.561521549276986, 7.284505096160084, 7.031618566660312, 6.873349902536235, 6.964267628521567]),
array([6.902492751232037]),
array([6.789139088692871, 5.61586645127625]),
array([5.510206011578141, 6.8662750225767555, 6.194984572986592, 7.120541471900228, 6.348939084519456, 4.466166882179358, 4.650873486954357, 4.904523103353717, 4.091493310775373, 4.246782004261961, 4.923858563031331, 3.6598480202171864, 4.216813598606841, 3.8787524959646356, 5.210403621401265, 5.053123030908158, 2.754435289929924, 2.5173347889091673, 2.34371663684326]),
array([5.575981405410968, 6.532080135294067, 6.549303188932172]),
array([6.982253211569698]),
array([7.109992181489436]),
array([6.376419036285524, 5.364997184952484, 5.595469816271385, 6.580675878970065, 5.76427145136763, 7.023498325337287, 5.552025664714753, 6.460257477706639, 6.915293244232582, 3.6778157497825816, 4.706790277983927, 4.860847832333356, 3.644950555519917, 4.859342449678877, 4.123761065370129, 4.863528380199549, 4.205569876220923, 4.374920718030917, 4.594142198786264, 2.6486432300560576, 2.792973227983664, 2.141024771453974, 2.1027962826112407, 2.5561215803963426, 1.2192689448135734]),
array([6.258950743979493, 6.936148753343894, 6.51926546783563, 6.731833445781766, 6.086724778285772]),
array([6.335239068791449, 6.539098285292937, 5.403840606470247, 5.94560782408022, 5.5264681022547215, 5.897817631320108, 5.883754975428849, 6.060304874052044, 5.199986837141498, 5.1225791693036244]),
array([6.343306278046706, 6.460551527201492, 6.820676236936637]),
array([5.411208192853552, 6.677442261542878, 6.680482150081403, 6.466739449317207, 5.806619349114019, 6.69013159571511, 6.342724448923436, 5.934426339480749, 5.810431672234408, 3.920187291023408, 3.8223388239174145, 4.638794912853006, 4.063171025156231, 3.7737645542983684, 2.5392489280425363, 1.777336395590459, 1.3028139083841555, 0.04629613905979163, 0.0]),
array([6.231552944903527, 4.553316240410253, 5.205094172839725, 0.0]),
array([5.524141121001163, 5.824812462019429, 6.071804357971126, 5.5830395510904705, 6.698508536777455, 6.583318142378004, 3.7617607873058354, 4.447464117155839, 4.3835154691040055, 3.8594741348344854, 3.9818663370619074, 4.366922784618428, 4.785775536582649, 5.068683761873082, 5.09186173200702, 3.45885558026096, 2.0769368646030615, 1.4730999911117928, 1.5166076133377584]),
array([5.901900616224285, 5.9383250414808195]),
array([6.526315547567774, 6.53267303948149]),
array([6.2125423466775365, 6.258742522155272, 6.1980109881387015, 6.581754997870371, 6.175061986922791, 5.421870073858536, 3.993622710791475, 3.614991263954094, 3.5701297584352054, 0.0]),
array([5.763955662640651]),
array([3.8488070266608423, 4.12697502033669, 4.432885911459651, 5.247220888646781, 3.3663388993833325, 1.2746528815358262]),
array([5.60588212792537, 6.170918374943317, 5.826456576153668, 6.214855008582594]),
array([5.982368924265457, 5.880491885213995, 3.78357785102716, 3.865370047587704, 3.6940935302563203, 3.0851813900874, 0.008935024018956628, 0.0]),
array([5.769217962886581, 5.732110100929852, 4.603158196727304, 4.3024854130636, 4.5208400011971595, 4.770333610217504]),
array([5.721351662642831, 5.752634248072457, 3.852211010574825, 3.8467758167335493, 4.536346681223106, 3.6703767793105753, 4.072104574242885, 2.76422076118822, 1.2038797813939581, 1.592410368066812, 0.5319921582899426, 0.04218138810291201, 0.0]),
array([5.428150471352007, 5.577562929254746, 5.451454394477877, 4.480852260117006, 3.751417659855009, 4.44220991153206, 3.679641794511305, 4.032628018771417, 4.365848569500308, 2.6648182419869215, 3.5446023983991917, 3.2810396527486576, 2.4882531907275482, 0.9767888389126868, 0.06786430211075388, 0.0]),
array([3.6420818111737416, 4.4577300335880805, 3.659788575236683, 4.409975608275175, 3.708786260543989, 4.324045465679406, 4.926154647695197, 4.963871478301598, 4.980585126853133, 4.838260176305558, 3.9544615543543298, 5.121313767788284, 5.002066513973817, 3.9267645825542674, 3.3329327602097782, 2.7490144261738743, 2.877768613514199, 1.510104814007394, 1.3839151297492627, 0.6640166158917121, 0.6137840021256866, 0.1131733186598689, 0.07314171593540786, 0.056571234645174745, 0.0]),
array([4.239606394344266, 4.432045396295983, 5.027186276839169, 3.9019870656448568, 2.654486971377111, 3.1099056681826798, 2.647494068968303, 1.8064262315058013, 2.476303368503329, 0.9525019207000632, 1.0472924933435936, 0.087711843731229, 0.0]),
array([4.750193114514868, 4.736931462393894]),
array([4.231779940040222, 3.8264845736106783, 4.739410734479973, 2.8464254814556447]),
array([4.816448419925018, 1.68085751091953, 0.0]),
array([2.657546566165937, 1.958176245769009, 2.480955930315741, 0.0]),
array([3.853190140817956, 4.526109719588449, 3.619450907186616, 4.229755366407074, 3.9604492648940295, 3.700690403398557, 3.552770530300055, 2.2279759971314035, 2.524024222322975, 1.989677254049359, 1.502819212320075, 1.1343362088840236]),
array([4.22112960624693, 4.273921722312831, 4.101812876067875, 4.051885580576146]),
array([3.6650962720635643, 4.3123109975185265, 3.19385297316192, 2.819635321391531]),
array([3.9594701525204985, 2.7809236370564663]),
array([4.283193818153238, 3.469374901107853, 3.3916788322904203]),
array([0.04940902586355511, 0.030573372191262813, 0.0]),
array([3.9482511093082184, 4.160222560292798, 3.5589860758696084, 0.0]),
array([4.243045581317044, 3.7481810003926253]),
array([3.9518697057336443, 4.060943123184487]),
array([4.02680998682257, 3.9451028099938363, 2.894167452511687, 2.1828169483036635]),
array([0.44997774210897973, 0.0]),
array([3.954373199486547, 3.659940734957882, 3.932197468772967, 3.080452638213778]),
array([3.71109204363544, 3.7197254663774046, 3.0497305206622247, 2.7432819760124403, 2.5431888863547134, 1.5584492484669854]),
array([3.105120136718169, 2.2160996247353033, 0.792267481611957, 0.8483222238655117, 1.7763018131373631, 1.7546848171642764, 0.15410837388805876, 0.1199235419422128, 0.09447050251846872, 0.0]),
array([2.9140730084705604, 2.9072853605725837, 2.9020263405372884, 1.8388848577429864, 1.5239145190046257, 0.9316271907153291, 1.4682167474177157, 0.7350069680460186, 0.43397396833648944, 0.09720373086449272, 0.0]),
array([3.3177131523731656, 1.2599593197817367, 0.11750189791461683, 0.0]),
array([2.8019973612829876, 1.1489458078698518, 0.08442669673851473, 0.08789527732384528, 0.050989603830179864, 0.0]),
array([3.0252492338924055, 2.6965972022097335, 1.8013478380651349, 0.3859689409086241, 0.6722918373609853, 0.07028954599773346, 0.0]),
array([2.297781677187745, 2.4167252424800396, 2.4401963477852835, 1.5854756251914015, 0.0]),
array([0.916237831565994]),
array([2.1321095205605882, 0.042276119655194955, 0.0]),
array([2.0601112100952506, 1.1285536198690962, 0.0]),
array([2.263582308847761, 2.2509826302252636, 1.5793183423610906, 0.8128645536594896, 1.5848573195063498, 1.251253624956966, 0.07846652243991326, 0.0]),
array([1.8727219404026034, 1.4689562233439268]),
array([1.6773013563073091, 0.044552737446933166, 0.0]),
array([1.2572945404023792, 1.334190965954299]),
array([1.7956470561213225, 0.6335628283290746, 0.08298495328429706, 0.0]),
array([0.008250582303288637, 0.0]),
array([0.9868580004866012, 0.0]),
array([0.4007489502440279, 0.7143036190223877, 0.0]),
array([1.1460806331974487, 0.0]),
array([0.8536748806944614, 0.18060180114141677, 0.0]),
array([0.8486069315262529, 0.44007470836769114, 0.35726095538062735, 0.11650309470686293, 0.0]),
array([0.15456175827742458, 0.6780916199410738, 0.0]),
array([0.3523433699965183, 0.05031307799406205, 0.07788376516209067, 0.04910305261988632, 0.0]),
array([0.059482581511536084, 0.0]),
array([0.07161745540272019, 0.0]),
array([0.26934756570965324, 0.0227518338080131, 0.11371021120813524, 0.0]),
array([0.0609471240949145, 0.0]),
array([0.2788011613807303, 0.0]),
array([0.12570742775517102, 0.04206823720531169, 0.0])
]
d = [data_1]
names = ["50"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T46', 'T48', 'T50', 'T53', 'T54', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T81', 'T83', 'T84', 'T86', 'T87', 'T88', 'T89', 'T91', 'T93', 'T94', 'T95', 'T96', 'T98', 'T100', 'T101', 'T102', 'T105', 'T106', 'T107', 'T108', 'T109', 'T110', 'T111', 'T112', 'T116', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'T130', 'T131', 'T133', 'T134', 'T136', 'T138', 'T139', 'T140', 'T142', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T157', 'T158', 'T159', 'T161', 'T162', 'T163', 'T167', 'T168', 'T170', 'T172', 'T173', 'T174', 'T176', 'T179', 'T180', 'T182', 'T183', 'T186', 'T187', 'T189', 'T190', 'T191', 'T192', 'T197', 'T202', 'T203', 'T204', 'T205', 'T206', 'T208']
def get_taxa_names(): return taxa_names