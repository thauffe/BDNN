#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.55002663855193, 31.58099563524633, 32.753368725378586, 34.10235881806799, 31.936326686388988, 31.730672146418453, 32.68969413756384, 30.51019243186724, 30.792826115534478, 30.846289718305883, 31.71313850544054, 33.3407767024105, 30.487419542922446, 31.087755047240254, 34.80999459364789, 33.04544429461264, 31.448463830210805, 31.636164368454246, 31.650935003380468, 31.01581319859689, 34.15336329318454, 34.46771969903739, 32.10469564143626, 30.30156638976104, 30.52263610481184, 32.53183691563246, 32.88873740282699]),
array([32.86054031487575, 34.08218765460701, 31.21627868157039, 31.290393570299706, 32.63001637334372]),
array([34.49474220387506, 32.47907783155044, 30.99329941211566, 30.922740561707187, 30.852525086449408]),
array([34.29150142003979, 34.35191322580862, 32.89696240652781, 33.17336302547268, 33.503576013782606]),
array([33.50032491189735, 28.358921844165113, 32.44021936106795, 31.757831263393975, 32.16914166791149, 30.911752351081038, 30.61055156490007, 28.23046032041147, 31.6254286811554, 27.140965496620794, 26.637086636261696, 27.24145582005565, 27.81264422063297, 25.550414739483205]),
array([33.17329148386496]),
array([30.067490017890357, 28.86913449317942, 31.506791569528012, 31.132106231576024, 28.17091884089593, 31.6155203545322, 31.559250288652414, 29.693407736752814, 28.481556647734312, 29.45276945977209, 32.37715446102765, 32.40392723176979, 28.98491405989722, 32.42539713218978, 29.462291673705113, 32.325772353112015, 29.87792744998519, 23.733330567004568, 25.53503218078247, 26.796359731387735, 27.021094113375653, 24.91533887602128, 23.59264068010715, 27.247395100818938, 25.05768888059815, 24.708464492437027, 22.694703026431156, 22.755255068497238, 20.023178205670686, 20.273526693808495, 17.571552080649845, 18.90857211933391]),
array([31.9970430404976, 32.560480055809954, 32.66578253301943]),
array([29.36563978873692, 28.114455237101573, 30.904654268812976, 29.915980089642183, 27.60817751785985, 27.933443641579103]),
array([32.339400931133824, 31.485345330260287, 32.39234978696174, 31.949711274532945]),
array([29.226168568385447, 30.650161423699345, 29.359584196685635, 29.38904882889631, 30.343355175010604, 30.59868868167688, 30.641554084257997, 30.81002493998832, 25.656391224259444, 28.002187715267887, 23.10089036292131, 24.221461839799115, 26.91806399670249, 27.942597813769545, 23.346968805527492, 25.699728411671348, 21.755918895233]),
array([28.164727402370055, 28.226779384034202, 28.667296562287145, 27.889565552076167, 25.37423395710279, 27.163034203855997, 25.953790550736194]),
array([29.66547957407307]),
array([26.29373427053948, 27.769097023742912, 27.425160934327362]),
array([29.36281388748321, 28.420777091669482, 29.098092891899004, 25.8317277595768, 24.955356601951067, 26.40829329315822]),
array([28.659177637465007]),
array([28.99572397120273, 28.631225910862753, 29.287816220050562, 29.223009985404527, 28.465522827898845, 26.09498930626725, 26.80009059276637, 26.685106931922327]),
array([25.8135599761047, 23.159287224843798, 26.80886844096173, 25.003141508196745, 22.461869446030654, 19.4296956142303]),
array([28.666621325967526, 28.19512226286122, 28.95584157673787, 28.830078019295314, 27.860497041361015]),
array([28.579311169366406, 28.79010797535283, 24.990826321812456, 27.116397436535216, 24.133631316439086, 21.1513299912479]),
array([27.596001476122133]),
array([27.96570015095893]),
array([26.26767011489252, 27.575379400240728, 24.745849801078005, 27.04290902128784, 25.183040567963516, 25.217453914475783]),
array([24.992372810734143, 25.608058151695587, 26.130796532245263, 24.395768737668025, 22.68648851460962, 22.664661253511788, 21.105025829968653, 17.000220462285757, 17.020738382931615, 18.622584214302215, 16.155743524546317, 17.467425812746562, 20.255899643483037, 17.49275440287752, 20.093832341197288]),
array([26.86841489121647, 25.09870833128319, 26.735484996283475, 26.968509222515124, 26.54090518482441]),
array([26.143293865448477, 24.86736798646499]),
array([24.25088696162094, 24.400706613696002, 26.25215343567042, 23.522697603262014, 23.83267934518739, 25.807816059774428, 26.351426769148087, 26.130402979550045, 24.248807035790026, 22.626234229939456, 20.235744730007944, 17.795490563989347, 18.201922907500457, 18.575190122341187, 18.613078837244167, 19.04860890404035, 20.406068587720664, 19.491790301196293, 19.657875572053428, 18.736848722103606, 17.590336829274754, 17.422186586438347, 17.742256482570955, 17.798093236847798]),
array([25.061009004632147, 23.77440140675521, 24.79593101857952, 24.34639981884944, 25.84603572403442, 24.15717077582461, 24.75557227629212, 26.292242561397696, 20.266779025807402]),
array([24.703448179631057, 25.35369256649642, 25.771478967467278, 23.770185887072188, 21.45588803598965, 20.997015246275787, 22.214665518048463, 16.936228999998036, 18.902970786075194, 17.463004751918074, 16.77488785773636, 16.762415851577835, 16.103476882701965]),
array([24.963535856396167, 23.851638484262413, 26.12325998884093]),
array([23.127112515775966, 23.861203455899116, 26.16437204629986, 25.766278196288933, 25.243592159914584, 24.44992841994563]),
array([25.178720299932984, 25.64102700619261, 24.512442855433985, 20.98176473492557]),
array([25.24515432456425]),
array([25.548578074774586, 25.15247479851213]),
array([24.19637586337148, 23.363560139476377]),
array([24.032523735819726]),
array([23.820956967691533]),
array([20.943386705688987, 20.040221828501448, 15.38851341905812, 12.540523773231504, 9.489335798157262]),
array([24.22882770531575, 23.121792120869586]),
array([24.051856954586654, 20.168869159386308]),
array([23.71893167357105, 22.310705752616528]),
array([23.295716657536822, 21.076024821803642, 21.616614184587615, 20.18051215939262, 18.84285967738573, 18.666214229085078, 18.859090152616467, 19.014448772238815, 19.93174121075496, 19.445531195965664, 18.78796042078623, 18.512018205539448, 18.888583445921547, 18.86153009925198]),
array([21.452797735004488, 18.369327512827084, 18.705475256577486, 17.705203683423328, 17.665328093079353]),
array([23.0916534968049]),
array([20.65005327012132]),
array([23.191400661476106, 22.851234619241026, 16.670601154023633, 18.73027223482058, 18.827941669919717, 17.6888280510349, 18.528877387254912, 18.649635755639007, 19.289070148595737, 19.44257251684621, 16.125224339250742, 16.406216918738608, 17.207321885943358, 18.09881961113516, 18.113850230353805, 16.453245640528557, 20.12963124683939, 18.938362891473176, 17.179583026968483, 19.112919556356356, 18.139725590033244]),
array([22.87489278063813, 21.867857982396277, 21.185392174936645, 22.61846519168906, 20.621144388658962]),
array([21.528759828051435, 19.06864021844259]),
array([18.990481231464447, 19.647845831805114, 17.8152627883412, 17.512538388421206, 18.936009112206545, 19.234339278944873, 19.195324219830763]),
array([22.209544119419466, 21.721858951860746, 18.40590088064258, 16.873870129175806, 17.94008152793395, 18.657599804067942, 14.9305665251842, 14.709617989959488, 11.261179394971654, 11.169803911436873, 9.339680256904867, 5.723381917166035, 6.766834424313507, 5.8201474832016755, 6.667509509094993, 4.444272191403905, 4.331850136216346, 4.522088877858933, 3.632238485416197, 3.763773873231913, 2.8077646154577214]),
array([22.355766893303585, 21.842109504296033, 19.60347178583768, 18.794216486614065, 19.029998981752247, 16.496776988118445, 19.104009324104183, 18.423364538069546, 19.19272898909619, 19.85726773486035, 19.247203045332643, 17.26811415075157, 16.798630646916482, 16.932499673861976, 19.692817538075147]),
array([21.685034560572646]),
array([16.56085132259459, 20.11061847881047, 16.763168625482827, 17.47980654094365, 16.769898338367533, 17.312036220775784, 16.949101251185528, 14.831422352146848]),
array([18.542964400211826, 20.260615363401353, 17.673646509913077, 17.24267110776591, 18.21882241876215, 19.307765768169567, 19.116358909472662, 18.64717447885313, 15.782761791639137, 15.535082586813186, 13.06851763670663, 11.685160296315853, 8.831286935171532, 10.446894397444726, 11.615654626139056]),
array([21.198715806579713, 20.78922863722454, 20.986865943015637, 21.31459154123488, 20.064265754737942, 16.020621029260884, 16.27244761330329, 16.425020930369588, 17.60587532221849, 17.820252101328364, 18.316000066425325, 16.942591909932464, 19.94891344077207, 15.851797551674096, 15.026381741395458]),
array([18.770299252610283, 18.50644502244929, 20.01578715258457, 18.56289536784282, 19.019255587508567, 18.34389581541516, 19.220620728159265]),
array([19.18364390654836, 18.560273982430868, 17.66937649465515, 19.633136088519436, 18.813932655429088, 17.900697452095006]),
array([20.284944346044792]),
array([19.723237533871842, 19.205737355817035]),
array([19.708093859392697, 18.696490529782324, 20.082693948150812, 8.556407445166812, 10.840897561155519, 6.940746659337129, 6.142570207293467, 6.44053577981762, 5.681136010889322, 0.4589071475628061, 0.0]),
array([18.354699795147976, 18.893096468481605, 18.377741918197565, 18.580789360659193, 19.808667035363555, 18.668475092343193]),
array([19.39497609841904, 19.32080664965174, 20.158781124944845]),
array([19.763089782834204, 19.547388716316153, 17.655126318418343, 19.788313134882127, 19.864852576889916, 18.90430355507395, 18.024501403977627, 18.091957118666887, 15.938704077981182, 13.571174126731835, 13.30410597649422, 12.550254275596926, 13.670310864841994, 12.642615151629396, 7.446773674174623, 8.062556280762152, 11.142923423907042, 9.654283590954567, 8.232348768940481, 5.92311660964539, 6.0956876971818295, 5.961575417365713, 6.5967316997418886, 7.188068717752811, 5.551786919731365, 5.945741422105887, 5.428069810959068, 6.751175510261159, 5.33293226651704, 3.644858179568497, 4.070841008247758, 4.199316576895798, 1.8585413041000027, 2.190174546713835, 1.0946261380539308, 0.9981948865512945, 1.3791577120146787, 1.4127540583266593, 1.7937258866023218, 1.3227631100041994, 1.2116897289627087, 0.16894341705636473, 0.13283468957501876, 0.4662026323484539, 0.3260289132491515, 0.20719422469800786, 0.7575462225199734, 0.20342087187950186, 0.2860674893664208, 0.08140695481360384, 0.0]),
array([19.496486899103047, 16.032256289905284, 18.368436198449476, 17.053278955944783, 17.124098250398145, 19.988297648059064]),
array([18.409369104253756, 16.814565964715882, 14.689260906824705]),
array([19.906831361308342]),
array([19.537072528670112, 19.156225236720406, 19.60694331708355, 19.535902255650697, 18.61916152436881, 19.4967689460365, 19.14993098158493, 19.03985961604636]),
array([19.446745960337427, 17.789056971137494, 19.258156629494962, 17.508410074201997, 17.38774211642893, 18.486471445262715, 19.318993644136267, 18.482604157701857]),
array([18.525174692817632, 17.991261435018732, 12.275021176520548, 11.059580411163402, 10.690132993509122, 6.1572261564034525]),
array([18.73287502050235, 18.23048266466321, 17.677617421799365, 18.02826592427805, 16.488069949064162, 18.064480605429615, 17.257613204213243]),
array([17.901858345945552, 18.00045934692363, 18.022615599719153, 18.24297686610806]),
array([17.88725707302118, 17.57921195057977]),
array([17.838276367647463, 16.8818588569711, 16.127825343518854, 14.91394729654325, 14.427637896154637, 12.816939027557925]),
array([16.41296164984634, 14.39669661200315, 12.346823170488374, 12.280299464664834, 13.666144329018309, 13.708556505898596, 13.741356603381726, 7.643830829153726, 8.21431149312659, 5.666460717828333]),
array([16.09605565706248, 17.434232841047645, 17.12315876140027, 16.985747988847503, 15.500610552279227, 14.589396190404537, 13.972740995777396, 15.069070936463191, 13.454142407440369, 12.335253439124502, 12.24168342770363, 13.73685056420276, 12.96034406892468, 12.064633166976488, 13.336906994725185, 13.663170779807785, 7.943505926111772, 10.45946642397921, 7.850966894861273, 9.299720615505594, 10.05078016015049, 7.55506886480415, 9.397687780152896, 11.14865709277879, 8.202519085535968, 7.613343739860534, 7.468173438543322, 6.198735981411794, 7.14223996155676, 6.2215414963445905, 7.1489083641547415, 6.732013370118148, 6.815318524224914, 6.637921210546035, 7.015605201942423, 5.4761954239804504]),
array([16.16570215850074, 17.45538337081123, 16.34809586661189, 16.253327034946082, 17.84117081990003, 17.50093084486944]),
array([17.74374076009076]),
array([16.129926618351238, 17.6631933897538, 15.672945940710463, 13.134689856820287, 12.150043109361427, 12.466414126545706, 13.076612081125676, 13.762100771870598, 7.891756311688056, 9.140904909151573, 8.555308362951637, 7.608156034567286, 7.432554160409394, 8.528723111588723, 7.207452082255328, 7.053746372146807]),
array([16.69666284313789, 15.111184153507159]),
array([17.608990793567223, 17.69526722466472, 17.307684925299647, 17.693085143399045]),
array([16.402548224954455, 17.58683334831667, 13.165136566222095, 13.771346718201215, 12.586687275712313, 12.855198984077568, 12.749103845642106, 10.870471631739601, 9.08516876899499, 11.285801535692531, 9.065263253378399, 6.7444923192544985, 4.090478836203758, 4.163950226053127, 5.0605793111160216, 5.02077031384824, 4.0296719303187505, 2.466447671030259, 1.774774514474122, 0.0]),
array([16.677676296026196, 16.76863670665691, 17.266888233979163, 17.495997829475975]),
array([16.07947967273976, 16.19615337278658, 16.333979954134833, 12.961621863742678, 13.69775537119269, 13.209924364618972, 13.784296502969179]),
array([16.572234120468753, 14.510059531731999, 13.136628000098376, 12.789956030703108, 12.625340237936689, 12.204275677887908, 11.520289657563795]),
array([16.99405563922991]),
array([16.4858146569634, 15.964777263692161, 12.610103149947866, 13.453373087097653]),
array([16.38127529793725, 15.594924081841686]),
array([16.096323248748043, 16.20428036390045, 16.245667073293244, 16.215099096759907, 14.930728876342258]),
array([16.11359211105586, 15.042610974595872]),
array([15.06924668676062, 15.65533203635698, 13.851776517837553, 13.26981243948019, 13.694640843128402, 11.87136765489982, 8.556214023757903, 9.23571232471201, 9.42888245281511, 8.177866166547716, 6.926603871696323, 5.9856660853867965, 7.1594362176583815, 6.967052600647157, 5.775905335834505, 4.184235271966156, 4.335593659455807, 1.4355987890989994, 0.0]),
array([15.100762688262963]),
array([14.761715643235245, 12.607605163644573, 11.958695407067449, 13.65974020708745, 12.762106446156062, 11.830157483986092]),
array([14.403139989591367, 13.561985994076155, 13.511722061870302, 12.322875125098836, 12.460204532078672, 12.750969704693617]),
array([12.033391759181956, 12.736091405039717, 11.588398852522303]),
array([15.376713725773419]),
array([14.788696649564914, 14.166497321511528, 15.2596064751823, 13.41054071252934, 11.932417941699853, 11.749667869344583, 13.07581361354365, 9.384925029574172, 10.348370015022496, 7.2415046923308015, 4.2538587767005565, 4.553293029424549, 5.302171648164308, 4.0395892318304405, 3.755787351004191, 4.971006471774894, 2.5261279677193302, 2.550189299779527, 0.0]),
array([15.01489452345379, 14.292668156445876, 13.10473484280178, 13.69736630813916, 13.03248737100662]),
array([14.844906343015335]),
array([12.2333432056078, 10.213980601394747, 11.033398590448119, 8.732297970059005, 10.589467105099684, 9.65601600497153, 9.494272823441609, 7.472235673420195, 6.199778508630318, 3.8531368296238044, 3.35239158147491, 3.3232853877479447, 0.7193425493089637, 0.591584222038119, 0.0]),
array([13.982739207780552, 12.74209679164796, 12.028515643136727, 11.751748528529161, 12.015571076351375, 12.509611780498245, 12.480850687868365, 12.754913287508309, 10.434859517313255, 10.869027907341769, 8.374655272883524, 10.981586825861037, 8.618224880574598, 5.452577280171695, 7.229902292626243, 6.414708501043041, 4.918227546984985, 1.4615387787856928, 1.2240220544031661, 0.705415719872148, 0.457350098389764, 0.0]),
array([14.165512066461034, 14.139896988319022, 13.542633538234565, 13.219925995108664, 13.793734701085166]),
array([13.865172510894071]),
array([11.814866872760469]),
array([12.231296115918738, 12.910610334082614]),
array([12.246194640786694, 12.64104588076257, 13.636726426372888, 13.72376419668997, 12.626066230108933, 13.75995683573898, 11.60656295017469, 8.32032218039134, 9.248890695405162, 10.90215625249437, 9.823335898946816, 10.703656622409662, 11.059857102293758, 8.40239991564172, 11.206800959021153, 10.942946881970531]),
array([13.221494202186628, 13.342450089436934, 12.755076850372411, 12.548592251208715, 9.773671725101115]),
array([12.366338997197042, 7.815831718270564]),
array([11.908100416794424, 13.014987445048734, 12.355606213104764, 13.194337335610056, 9.991236860727168]),
array([12.750266945886523, 12.070651960962003, 12.586734053552714, 11.559427344773948]),
array([11.770187745924817, 10.345886686575104, 10.781828522554108, 11.298549840326992]),
array([12.709644989042488]),
array([12.277382932132666, 12.099793253191027, 11.726780249427929, 11.852584066615506, 12.44022127077352, 7.3151325061156, 10.91817707871628, 11.243772397124618, 11.052502274719142, 9.4895578014574, 9.907791052016519, 6.304505545769393, 7.093371596602946, 6.086332010083686, 5.449084621577816, 6.496550059200846, 3.848809430779, 4.326342418844221, 5.179633734728482, 3.3665756449505055, 3.2100055851089175, 1.1209839059115854, 1.068879333051958, 1.7583060403589859, 0.8170306561433681, 0.543397640375654, 0.4357253204122797, 0.7783322844098204, 0.0]),
array([9.646868607204759, 10.871766198566046]),
array([10.29163121616601, 10.26629421201067]),
array([9.901027495893766, 10.24936012280381, 7.533052404059295, 9.033570661474071, 6.730839436117361, 5.999289286796816, 6.339374952416305, 5.1456120355096635, 2.6759358315788226, 3.07976738318567, 1.1773652190217532, 0.7513055661054436, 0.5811790976030424, 0.0]),
array([8.195389193208351, 10.66359948411543, 9.128496407916847, 10.239727015968805, 7.302282964238076, 10.930666545098356, 8.059293972640138, 6.528428020042174, 7.240815276064559, 6.4117724798715106, 5.733245504559964, 6.03882593840749, 6.012743969875928, 5.229198538204964, 5.120298588160221]),
array([9.49984616513554]),
array([9.735488190197266, 9.06770054165105]),
array([9.647816366038786]),
array([9.525475331521813, 9.031050673831263]),
array([9.864353922630713, 9.838508354895271, 9.638910530084576]),
array([7.8511297558526865, 9.287411557245756, 7.03251109895816, 5.450516971103529, 6.691202752397263, 5.080854028609921]),
array([9.86754038364464, 7.7767418240536985, 6.264076363729975, 5.891998798202408, 5.90275909462216, 4.1224066027170645, 5.104046300436348, 4.2751695156948655]),
array([8.874132669043718]),
array([9.198195121325364]),
array([8.991719428838536, 7.186760628529642, 6.115026934424925, 6.067730053825516, 5.282923245705396, 4.942653077140718, 4.376501296383877, 3.6115311366846248, 3.1123327152669544, 1.4554056780834548, 1.3881445518697628, 0.0]),
array([8.42428889188274, 7.456606241139186, 5.792392790801207, 5.45072425167329, 6.577163940628752, 7.240828704055081, 7.177178972007516, 6.390473239035522, 6.889518048991335, 4.044145915251981, 5.089632767892749, 3.0770808082392884, 2.3902728246345437, 1.5177048734675005, 0.2847907725895374, 0.22498883849943518, 0.2814834886666483, 0.04278064916657162, 0.0]),
array([7.977964556036162, 8.173689490081259, 6.339210699838246, 6.957479089082202, 5.77757667751224, 6.1096117168219894, 6.446543683808031, 6.760452607037055, 6.130433692907105, 7.170648174888084, 6.436729736524594, 6.970382534592767, 4.751820893844628, 4.106858059404452, 4.575703321756914, 4.149338719961679, 5.071749407166219, 4.92777242075315, 2.7988592512753345, 3.4992435063446767, 3.16359859509869, 2.8157605523579554, 1.6685676020506381]),
array([8.056791807345524, 7.28843984657468, 5.830792174725229, 6.64163132119714, 6.650283899315986, 5.142952733287633, 4.336579772688049, 3.268823625571792]),
array([7.869968248594026, 7.3474002632465965, 5.508214328685602, 4.893833285278198]),
array([5.694408833342221, 6.390310681058834, 4.212740643067977, 1.516128468476182, 0.7187001875637576, 0.7531822804742656, 0.5484281124839259]),
array([5.896647955919051, 6.862485640858182, 5.715674021863495, 4.702715115995485, 4.574566621678831, 3.7423945674840735, 4.502758994869841]),
array([7.651201500328027, 6.389903767355441, 5.963939858679487, 6.145744160237875, 6.561005197938381, 6.519808753061474, 3.7926833698128632, 4.193961313533334]),
array([2.4009732682498175, 1.697372221682362, 0.39258785019403586, 0.0]),
array([6.081923690940228, 6.444355962297172, 6.072649702158323, 6.3396603048820594]),
array([5.405044930500891, 5.967789393375609, 5.455073353869707, 6.852049438348744, 5.666623998498243, 5.465456635115143, 6.839051016046371, 6.184605203203828, 5.781253298477097, 6.1979370600598624, 6.65023457656464, 6.004329479242334, 4.349813604892255, 4.879819506726206, 4.159901736213732]),
array([6.802304019918797, 6.610738176888873, 6.356185583878579, 5.8135906698528395, 6.166856005483734, 5.524703010309437]),
array([6.727515519355356]),
array([6.388791763558633, 5.4352683510588236]),
array([6.034274585074007, 5.502777035125513, 6.245298902841034, 6.233230075684966, 6.1287276404340645, 4.8957137856144115]),
array([6.196272058444561, 5.638645456217408, 6.063289783669391, 5.789657854458893, 4.448950689336826, 3.668944516989895, 3.4078572074591014, 1.5254261283163462, 1.3105131015832394, 0.819039662058019, 0.18076954516511412, 0.3259978356765054, 0.05091851808479586, 0.0]),
array([5.857059644850573, 6.065563127993393, 6.213360516942158, 3.9032228309007007, 4.488267351615682, 4.968611529188704, 3.107897404387984]),
array([5.850529369387591, 6.173956750125475, 5.200725079298246]),
array([5.502997035326675, 6.296217751840569, 3.7529998646992566, 2.4198492508989755, 0.05073681032337682, 0.0]),
array([5.441842402340022, 5.361911391689213, 4.8929200581915735, 3.8163609863529633, 5.278641802656496, 3.9908455542682884, 4.234274762963237, 4.464742048121337, 4.400514156193241, 4.802392726257284, 3.113551310172598]),
array([5.333748350471938, 5.368143515844567, 5.46769565653042]),
array([3.953244984244714, 4.0232354982473755, 4.08559437730634]),
array([3.6307787829534086, 5.274674891253154, 3.5242488421644755]),
array([5.465272010943696, 5.256096214376515, 4.85890225452939, 4.046879279500186, 4.354557313163791, 5.00726289173862, 4.461562685845942, 4.066494015152326, 2.7863721020785617, 2.970653039315919, 2.0482436757190996, 2.1755531730888844, 0.9328024043669447, 1.0163767668259653, 1.026312991937836, 1.0858194645440222, 0.5705908446534915, 0.1830137080526888, 0.2498355858991912, 0.6625770522756084, 0.2833221624047729, 0.0]),
array([4.90367096718126, 3.9624514306999292, 5.00219728139862]),
array([1.2706418793010297]),
array([4.138285447268174, 4.7710952459853635, 4.758803575964034, 3.876088492731326, 2.73482062747554, 2.9979028769460667, 0.9896136041966397, 1.1633277591227342, 1.2573657225623607, 1.3102751722612205]),
array([4.544585540786828]),
array([4.120573478938427, 3.663868076823268, 3.4601015595967803, 3.184976340180801, 0.6081458432640635, 0.08378861667836085, 0.0]),
array([3.7075005454651278, 3.8011014767338533, 3.607579160021106, 3.4653156088550663, 3.251576443036165, 3.4271855467917165, 2.5208693740233565, 1.786412394455929, 1.4065159938982332]),
array([4.121212632387253]),
array([3.0410346158063204]),
array([2.3083972031146422, 2.0485783378245137, 0.7533788438850336, 0.2753264125999112, 0.0]),
array([2.327283491693018]),
array([1.5157428220873597, 0.30471161141068404, 0.6540540503757966, 0.28325635276063293, 0.5301233399106776, 0.0]),
array([0.5910183890443338, 0.5396359089997792, 0.0]),
array([3.319433693921208]),
array([2.396751031845955, 0.13195520307224573, 0.2802863235443198, 0.0]),
array([3.156050574884178, 0.9634816080679842, 0.9804425498091641, 1.2718231266118432, 0.0]),
array([2.657210233330879, 1.3675715122008378, 1.2257666430024874, 0.33877981844108157, 0.6012154139063819, 0.44821566078198943, 0.48332148493088717, 0.0]),
array([2.4139021863258954, 0.9923328264037012, 0.0]),
array([2.5939179669861736, 1.6577231667470507, 1.6663513371987786, 1.1236907065359696, 1.4045888907638853, 0.15858036074025317, 0.29452320138330623, 0.06312116899835238, 0.0]),
array([2.291062469926567, 1.1047388363595565, 0.9231817281500121, 0.17569580701758558, 0.20742874955462454, 0.0]),
array([2.9728848432285693, 2.6296557811167154, 2.4602420942253973, 1.011281427758866, 0.7679160856107713, 0.5016114684932482, 0.7281898917792476, 0.0]),
array([2.7881168200683257, 3.122834886592114]),
array([2.525232460934033, 1.941378584822793, 1.29422766921497, 0.6568396393060981, 0.4972642019177086, 0.03618987720326239, 0.0]),
array([1.9928160662608412, 1.2472881127738211]),
array([1.1833288886093603, 1.2600712696307421, 1.1712927337582493, 0.16065791466187662, 0.0]),
array([0.4027193869423757, 0.0]),
array([2.5615261976254615]),
array([0.4395902710930917, 0.0]),
array([2.0385017461591475, 0.17588971117014895, 0.2056189793208536, 0.38849513089958687, 0.024464233882737885, 0.0]),
array([1.4157368033642324, 1.6021446861540372]),
array([1.4822885491666504, 1.085101462485596, 1.2244162221340067, 1.6835320180910303, 0.6735697945414872, 0.1641292319916715, 0.10721426266806486, 0.0]),
array([2.0529190804373]),
array([2.0388265100239757, 0.8938988365623182, 1.573779313422243, 1.0677121314872147]),
array([1.2437068498060535, 1.2102053901041268, 0.24788074763555523, 0.0]),
array([2.212786659592061, 0.4358919260070472]),
array([1.676983746344066, 1.3511997111414302, 1.5634605340287282, 1.5994389425826436, 0.2262257609932823, 0.14970617283002075, 0.4047924398988893, 0.7428166171047046, 0.0]),
array([1.0549628849696804, 1.1140294560360946, 1.7073750623155246, 0.8853415633055571, 1.572092783002523, 0.004743146437781204, 0.0]),
array([1.2492282469563314]),
array([0.7427230862106451, 0.41807873377832955, 0.5100015216184697]),
array([1.1281177921954697, 1.3343295766033414, 0.5835950956509418, 0.0]),
array([1.5557438155583647]),
array([1.4233970733527457, 0.23869450570573691, 0.0]),
array([1.463385179738508]),
array([0.11954227228258507, 0.0]),
array([1.5681324250529816, 0.8087784512569469, 1.7099169466263098, 1.1046929035589164, 1.5401004910007898, 0.10138521893487108, 0.0]),
array([1.4994104613635544]),
array([0.22701763343527503]),
array([1.226411959783994, 1.1991281719316271, 1.3542105684499273, 0.2925087364057449, 0.6710486972902411, 0.6214642252238152]),
array([0.3634225194086832, 0.27207787327493016, 0.12791077436248788, 0.4226814789579743, 0.0]),
array([0.875311515406617, 0.0]),
array([0.6799194599680891, 0.40546419163305863, 0.5242626528449696, 0.0]),
array([0.8333427275994381, 0.364971074347188, 0.19400737294670278, 0.7512271572362305, 0.0]),
array([0.9061496524002477, 0.7353112962983116, 0.44419171481474135, 0.24342594923009453, 0.0]),
array([0.206641962375176, 0.30671759825227557, 0.0]),
array([0.3651086150218785, 0.0]),
array([0.2687043150650279, 0.2989831048977076, 0.0]),
array([0.3116287098330873, 0.0900789230411805, 0.0])
]
d = [data_1]
names = ["63"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T43', 'T45', 'T46', 'T48', 'T49', 'T51', 'T53', 'T54', 'T55', 'T56', 'T58', 'T60', 'T61', 'T62', 'T63', 'T65', 'T66', 'T68', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T79', 'T80', 'T81', 'T83', 'T84', 'T85', 'T86', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T93', 'T94', 'T95', 'T96', 'T98', 'T100', 'T101', 'T105', 'T106', 'T107', 'T108', 'T109', 'T110', 'T111', 'T112', 'T113', 'T114', 'T115', 'T116', 'T118', 'T120', 'T121', 'T122', 'T124', 'T125', 'T126', 'T128', 'T129', 'T130', 'T131', 'T132', 'T133', 'T134', 'T136', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T146', 'T147', 'T148', 'T149', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T160', 'T163', 'T164', 'T166', 'T168', 'T169', 'T170', 'T171', 'T172', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T179', 'T180', 'T181', 'T184', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T195', 'T196', 'T197', 'T198', 'T199', 'T200', 'T201', 'T202', 'T203', 'T204', 'T206', 'T207', 'T209', 'T211', 'T212', 'T214', 'T215', 'T216', 'T217', 'T219', 'T221', 'T222', 'T224', 'T225', 'T227', 'T228', 'T230', 'T231', 'T232', 'T233', 'T234', 'T235', 'T236', 'T237', 'T239', 'T240', 'T241', 'T242', 'T243', 'T244', 'T246', 'T247', 'T248', 'T249', 'T250', 'T252', 'T253', 'T254', 'T255']
def get_taxa_names(): return taxa_names