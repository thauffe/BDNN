#!/usr/bin/env python
from numpy import *
data_1 = [
array([22.447354752041505]),
array([30.887382229223363, 30.21491929653252, 31.104484526665598, 28.402232176524457, 33.54258425625262, 25.197722612871537, 24.937471326147968, 25.133141103755136, 25.5727963644458, 25.374506648633524, 22.81647997159976, 20.98513594053625, 21.177231040890852, 22.221788412871497, 20.926340730851273, 19.070001731086332, 15.840761857763306, 15.7073223016059, 13.880979076065058, 12.220466093231336, 13.301408319854787, 12.012282955896296, 11.133116986233254, 6.143732686672166, 5.246957848153159, 3.882553901949585]),
array([31.281078942024777, 25.97595431267885, 26.804045348474432, 24.96747867052739, 27.284974310298434]),
array([31.47776503447841, 28.80084299430867, 29.649123507563594, 28.186362862322184, 26.466998878011335, 27.046134886568986, 22.192021215427303, 20.033990120872254]),
array([29.904697734643506, 30.38316479457907, 28.16650922630278, 29.21630207383268, 29.888711880783188, 29.11860407003374, 26.746083671411686, 23.543317718501488, 24.239592456961738, 27.318089446855495, 26.317786769278825, 26.148109325465057, 25.9605417804181, 24.142924051843266, 21.876051256836508, 22.347498458925706, 21.536093536598234, 22.299536114830378, 21.949668418649644, 21.352169988077463, 20.821094483648675, 20.675271372194565, 18.02866689733989]),
array([28.52869078508203, 29.057234711970274, 29.46393558434695, 28.3679375059708, 25.043541949291562, 26.589900127917137, 21.00029522273484, 22.19929451249027, 21.033305329458337, 19.020191885962433, 14.706424848782058, 15.596601128036472]),
array([28.19445895699349, 24.532461082351134, 25.32125700200169, 22.25556265918229, 20.815859315195766, 21.686855294953173, 22.35279129445002, 15.588976018645127, 15.213321202585016]),
array([28.41311468186443, 28.84403775421696, 28.771241610584326, 28.61971707225216, 28.83577272306465, 28.000904371437134, 25.75938413617668, 27.624700446857556]),
array([26.095319200157096, 23.073825918978542, 27.051645445077316, 24.51941711157341, 26.567398076691394, 23.19124703799769, 23.61859619488236, 26.583099388011433, 22.106011616238735, 22.24382036040257, 21.186770638002955, 22.595994256022106, 21.859099034105807, 21.26534215256648, 21.591354734651826, 21.867605661012256, 19.965632907641968, 15.022645178746405, 14.04852707982705, 13.897318673426915, 15.76149058362325, 14.693180371994302, 15.365372711440923, 14.371436823015419, 15.550733525727862, 13.052778429282235, 13.742037182539695, 10.081985993303256, 9.943872318206328, 8.418674805468925, 11.137748216818872, 6.925461494225313, 6.067727998740767, 4.358777814949375, 2.23302505681263, 1.3478645203834663, 0.8106496419657816, 1.5659819881880588, 0.0]),
array([22.392112683378443, 23.011773257381435]),
array([25.66650217306197, 25.418502166780264]),
array([24.273670869390795, 24.002253071230534, 21.494660086797374, 21.50084940298445, 22.752545312037594, 19.621178850300733]),
array([24.583858086613677]),
array([23.659076085041107, 20.682160327749465, 21.66009898446791, 22.72954636433276, 21.11959936911907, 21.72106746124736, 22.61729135470916, 17.72377829230846, 14.161470842507025, 15.102851459665523, 12.328643352811307]),
array([23.177714304040375, 23.43832957136455, 21.615026580143002, 21.797275861291954, 20.444952395994513, 22.236018875666367, 20.455047898536094, 20.65298828834646, 20.532172267178627, 22.131252679500847, 21.136604569012537, 21.676002438848656, 21.780911062868217, 21.29074314306712, 23.015735249146516, 21.12246418941679, 21.19650260209131, 20.471001218437976, 22.680551402826662, 21.640112165537257, 20.69575138034515, 20.972853900828895, 21.50481226337622, 21.181280048167476, 21.421621951519583, 20.826396191956377, 21.513398322737928, 22.870850588854502, 22.41523189953889, 18.657185027848573, 18.783134894190614, 19.252096468801987, 17.701037704604392, 16.923811971707643, 15.706240905754882, 15.441866033156938, 15.492814362688627, 15.005470229254094, 13.906542588109449, 14.185593306280493, 15.450916415807814, 14.66757106291379, 13.353172920594881, 11.9578355144622, 11.65072107018148, 13.380378472242263, 12.577439122849714, 7.824844747847332, 11.211802074682284, 7.431669427184327, 7.926166829816385, 6.923790592618424, 6.7540941625662665, 7.06892231919675]),
array([22.152047165308808, 20.501006530691637, 21.48775682146337, 22.534253072928024, 21.295657712407248, 21.44665837027246, 17.321914982312734, 17.69737088846724, 15.325358639509076]),
array([14.947900981495966, 14.755439759157873, 7.491111041736673, 3.807975666059363, 4.077018212861775, 1.8310506556919568, 1.8243573137025084, 0.0]),
array([21.849169223279457, 21.41844622125209, 22.486808893739873, 21.679944002115022, 20.92556926237557, 21.107850938123708, 22.071524509905146, 21.436934718705306, 20.93670311712882, 22.232690238102492, 21.256467248612093, 21.69712817597086, 20.727318865359493, 18.09187588941493, 19.433588216327824, 18.18301370631015, 20.166113446074004, 18.580502129435484, 15.154887611295415, 15.712259239705013, 14.923481457690709, 14.156405888188818, 15.194104328910857, 15.724132484524134, 13.948134234459967, 11.668455963672777, 13.051752012722247, 12.868998541548619, 8.600827479437005, 8.953285893009546, 10.189768672885766, 9.81766740458349, 8.278965612327315, 10.309618787853164, 9.724339702048821, 9.361728854448858, 6.10903816885272, 6.9840186354484, 6.837639815265111, 3.869191050978813, 4.244007031086372, 2.2261955694023894, 1.5546833131745739, 1.5151932878671401, 1.5326515726501109, 1.177007507486084]),
array([6.76939096593812, 0.0]),
array([20.908045860584835, 20.48018133119071, 21.661345803946713, 21.45031126835823, 21.445369998252534, 22.178342059835106, 21.777904006031633, 21.308161273534065, 20.49483209139646, 21.613877143503252, 21.48447595225683, 21.288939697702055, 21.1931890790312, 20.752805318956625, 22.226964872304794, 21.119948219819218, 20.9619424431937, 20.66659775442037, 21.276515112013556, 20.92330282358852, 22.09716302114407, 21.214886033579226, 20.919347453482224, 19.78010947544857, 19.604989791356484, 16.0956841038742, 18.682057741366712, 17.45559160265842, 19.471764922693033, 18.59431288862553, 14.179672605182915, 13.922332923035174, 15.43484982700789, 15.367428485890391, 15.702892273482833, 15.708495512706588, 15.253633856460063, 15.820301869271871, 14.11867631976557, 15.174600230270256, 14.393733847913206, 14.921528984304793, 15.123854757074724, 15.830273210395376, 14.99804303639564, 15.675570560099427, 14.116158729962484, 14.400648021372717, 15.90372161725638, 14.719261811505225, 13.984220350694455, 15.966704885230552, 13.824830174373124, 13.950607684900046, 13.4335588406456, 12.700493036985868]),
array([20.74554385215697, 21.381038723237115, 21.083902721995734, 21.102354859622114, 20.764028137813305]),
array([20.514161402457336, 16.753676897078446]),
array([20.446798696280347, 20.455472151096554, 17.60911849060799, 17.1460109570804, 16.929378845798468, 19.764235141864642, 19.311929883776877, 18.716637733253293, 16.141286519337825, 18.10019935866435, 17.561095956241424, 19.62207332012188, 19.047594664800254, 17.701239388669617, 15.282463494349713, 13.911310984710937, 15.186183421004882, 14.68121175867325, 15.546764473407768, 15.028752967145486, 14.72569232960922, 14.865450223029647, 15.355738871780037, 14.95466533460402, 14.3804853899884, 14.433428050570775, 15.841909941133888, 15.250159327685013, 14.497089536674885, 14.36958619469203, 13.872624014518793, 15.374808528466039, 12.14800841057564, 13.552571254758112, 11.841082248983541, 13.2936758608776, 13.547381335673636, 12.188964915712928, 13.531621541543945, 12.144235047042823, 12.138084890632596, 13.730745424588934, 12.35863859811246, 12.181092529698823, 9.596888263444676, 9.255625672588392, 8.075001884383525, 7.729008172349501, 8.087264343838637, 10.071072308407395, 8.312578559901821, 8.31418564160171, 8.181850173315105, 8.003971372835577, 11.387022144905002, 5.840922696959274, 6.887703388028369, 6.738025642172241, 6.449837305211156, 6.3426443130813634, 6.560936295535127, 5.773057979347584, 6.328511822493792, 6.269028153408111, 4.454702782360597, 4.786562491263075, 3.908483546440263, 4.410654685440022, 4.55800111535932, 3.891208842255744, 4.776674122522206, 4.585079650465004, 4.228402465952893]),
array([16.34735103505466, 19.557733504182515, 18.891138304268672, 16.717810899354653, 19.88629186239967, 16.377554910673698, 18.60399902709733, 16.523543207579685, 17.622126551844488, 16.193144523351194, 14.47568184569225, 14.145385433332148, 14.970007243657994, 15.009044666166231, 15.423277023270014, 14.33543032501509, 15.671702860958499, 15.419259749693236, 14.22719173365806, 15.64995646290426, 15.766603720366112, 15.78307144055071, 14.64502891242138, 15.931716389056383, 14.43078852592362]),
array([15.216577694656015, 14.173526283615484, 10.302572107965663]),
array([17.715421559007968, 16.571573551996252, 13.963012598372176, 15.291155065772836, 14.264403671803864, 12.625415765908869, 13.156968761927661]),
array([15.974067694518174, 14.497177239081845, 14.750878374605332]),
array([19.184762154886645, 18.32679808092592, 17.932837053387978, 18.615621725389833, 17.3499852539488, 15.226454970472956, 14.774095958662267, 15.827691437091632, 14.278714454802236, 14.267005517099008, 15.462429911512908, 15.402677930391274, 14.618905346682036, 15.001805080737753, 15.620042793276294, 14.99201020516946, 15.795745299310477, 15.765709829008022, 14.116813022761718, 15.929122995074485, 15.415545432134184, 15.251728114859722, 14.163081501516158, 12.15376858362173, 12.369332892632595, 13.733837182953266, 13.52314068037731, 11.239553271012863, 10.576009620424129, 10.74749551204481, 9.918414244021545, 10.33935185119276, 7.38165296881148, 10.791506537204292, 7.924000968132377, 10.847829060312387, 9.670694039858542, 5.6519507836700305, 5.471118259890594, 5.34197253350295, 4.758810886212609, 4.840185839153166, 4.12556207762462, 3.8684436564287905, 4.806973553226616, 4.2871095691890755, 4.7234487203533595, 3.326567526097, 3.4102104476884683]),
array([17.546175323647134, 16.150656364104595, 13.854842243243956, 14.291633375811783, 15.63369294977696, 13.517210019299577, 11.56943040776872, 9.434563128731558, 6.106848592605557, 6.0558899813570255, 4.533448957284844, 4.476796460312569, 4.77047355280433, 4.471554757826176, 4.0882739091879134, 3.5356521128522704, 3.1699896814030253, 1.8313913420366603, 2.5447797118986397, 2.4902455214732626, 1.8086905096623824, 0.879090172199026, 0.0]),
array([16.869012738474186, 13.894967925206117, 14.309317153825674, 14.469007238797449, 12.872857207825163, 12.500674358602076]),
array([17.12494258400138, 17.142709061106082, 17.214361427751196, 17.681198716777136, 16.257097095890895, 14.852878532309841, 15.45444640052197, 14.737412809252731, 15.44248317714872, 15.54123494994952, 14.145275245104774, 15.658558434713665, 15.01877392451119, 14.485873015915171, 15.135492060980305, 15.215273259398431, 15.266256698294153, 14.920207130817536, 15.864073608159247, 14.663629116726993, 14.376614782958335, 12.875096567531864, 11.83723303084305, 13.411615643288638, 12.726301758258456, 11.768248090463198, 12.334886967918337, 8.738865570018927, 7.480644996394059, 8.722008001309371, 10.75632972242802, 9.997669075752997, 10.348529485482263, 7.373517471119536, 11.44970542182217, 7.649268750505126, 7.638094169522109, 11.516922270676108, 11.287115353424696, 7.154533411450705, 7.028644119420762]),
array([16.74440193678976, 12.453004063255685, 11.06286352592818]),
array([16.258588372066036, 15.999230418976769, 16.797502508453427, 14.786193938928033, 15.058101684959986, 15.607638395016343, 15.023009933819871, 14.037955923900235, 13.97685766293309, 14.130712588033841, 14.313289223302865, 15.838137220973852, 14.993423426419358, 13.46592172280571, 12.917031516916632, 12.76299965066561, 13.392303204844042]),
array([16.42625092112149, 14.680508171660069, 15.73365844473943, 14.391065509429502, 10.579153453493454]),
array([14.650987998887963, 15.06827079204042, 15.430664956469103, 15.125913462646622, 12.367322275819697, 5.701395911773713, 7.066769347935988, 5.183758126701747, 2.0028526101262236, 2.164860551967681, 1.044268569460859, 0.6250026253428784, 0.0]),
array([15.015133045350892, 14.306262720016436, 14.263204586639494, 15.50746376600528, 13.918944679508918, 14.431414377556564, 14.518501835475773, 14.155996215035834, 15.17653106014228, 11.797805672108007, 13.271926813975805, 13.471448819481521, 13.565683797116012, 10.572283444336907, 10.512847037276575, 10.289831882407581]),
array([15.47613559734494, 15.529534325957616, 15.225216710834237, 15.763961964650099, 15.598196416738164, 15.239814779088018, 15.406636204020636, 15.561625389095372]),
array([14.61266000098919, 14.600091693047764, 13.845638558809165, 14.82923614274379, 14.827098900306936, 14.732281061956982, 13.919159117436708, 14.231871456907424, 13.896970345931804, 12.504508515373788, 13.81939045048344, 11.29856835944436, 10.544526049751209, 10.665176283837805, 11.038459796584222, 7.285318361896862, 9.764226238244033, 9.62741830707699, 5.504945555864227, 6.289067653382126, 6.426038150093912, 6.390356964174165, 5.242698101566045, 4.276317887367586, 5.027185010938459, 4.044972824920022, 4.5131641675652485, 2.217345224447535, 1.8331441617885642, 1.9731511567279165, 2.3432604424316477, 2.5144921548880483, 1.977672253749533, 2.081526985361008, 2.317776016331843, 1.7579360537058617, 0.8204745766533185, 0.9799159943599305, 1.736495993250867, 1.590146306667235, 1.6194022119752245, 0.7060212423513479, 0.747888079752852, 0.17809340447896616, 0.0]),
array([13.864169681373763, 14.030545916546606, 13.591381276924674]),
array([13.041756583842014, 12.411278977638384, 12.499445286557068, 11.282179383433608]),
array([11.9408439757026]),
array([13.97822039594411, 13.973874794729124]),
array([13.860380169367463, 13.934665567537762, 12.1586336709419, 13.085130732761751, 12.90108184395714, 13.788462705671524, 12.502144218677948, 13.726144420253476, 12.409373199322475, 12.55644877460729, 12.72373519862523, 12.670269275667188, 10.398722515258733, 9.002618287844257, 9.609543457652846, 8.415006225633721, 9.98792309749933, 10.452461763922942, 8.01976800437604, 9.903711874715615, 10.478959431387777, 10.510722478014607]),
array([11.005251086581081, 11.158476065095504, 1.8763938258815083, 2.057317890254246, 0.0]),
array([12.440451792554612, 13.467136884454396, 13.56694539169338, 13.414725675021316, 10.970694106780964, 11.133804903018573]),
array([12.510357463913437, 13.224196663392501, 11.815117140127605, 12.38646184528816, 11.285663367681835, 10.784245271774727, 10.836504132521418, 10.87765394881855]),
array([10.298267360054156]),
array([12.734450973055287, 12.936632277338571, 12.244867225013468]),
array([12.321467940415296]),
array([11.971024697086586, 12.076806042315962, 11.771914282432293, 12.562936988140748, 11.320352234599591, 10.655714926317348, 8.169896766316366, 8.168921946312867, 8.044864749235401, 11.567908691334097, 9.172125986521811, 11.314020298416548, 9.711601735680121, 8.344286655353677, 7.942757155804367, 7.414119367157364, 8.578889896015923, 10.5474619287499, 6.229389159783609, 6.704405912305923, 6.536936802069279, 5.518631430402644, 6.24585683810946, 5.230616502508794]),
array([11.869112308844814, 11.541600131209334]),
array([11.797020196232907, 11.819914622894078, 11.197889931293812, 7.886256644203776, 9.161052027245182, 9.931535819670454, 9.593492375698627, 11.26916885237648, 7.829413898211399, 7.422775647962374, 10.210407726593976, 10.817855219440267, 6.871306514076387, 4.157015703166996, 3.8564989299180668, 4.602267261937556, 5.281428842132993, 5.126763067761524, 4.047126239395336, 4.728218192892424, 4.178156843626127, 5.172525520289881, 4.2852498448409815, 4.076973868789365, 4.740441317668147, 4.621941513920005, 3.627818138841766, 3.650758025455538, 3.6779894830191466, 3.4158176689604067, 2.298171805186874, 2.040968846061749, 2.5052440235668154, 1.856849324131426, 2.5762823361115204, 2.193201087367619, 2.148089880186712, 2.399076830597243, 2.485229113769415, 2.5515046838135333, 2.095664360761753, 1.9631269715205248, 2.35909344785508, 2.4011459374251007, 1.9937534877209702, 2.5712939391675014, 1.0393168150048222, 1.7113644847212188, 1.6704618361338752, 1.65788772473039, 1.2950363332331924, 1.763680431074704, 1.7861749095142203, 1.2978912581500608, 1.1024674823938228, 1.7033460701846885, 0.9767561658679864, 1.4533734293963807, 1.0113008061648214, 1.0554904476421374, 1.794094050559278, 0.20396246377861926, 0.227676052686703, 0.23281899737608336, 0.45325876158582346, 0.6798634086669887, 0.0]),
array([11.796273563653322, 10.040281496123901]),
array([7.347827458064767, 8.188200430889292, 4.555126275124472, 4.775332948800135, 2.527809142588247, 0.0]),
array([10.53491369285808, 9.315768121301371, 7.6305243725346905, 7.427384000985946, 4.06356145098942, 4.403573285920626, 1.9299498703185725, 2.0660613747837875, 2.2102517283017282, 1.1568105947928333]),
array([11.796459348702538, 10.80835789298731, 8.894580156899847, 10.550968013584008]),
array([10.575137091385486, 11.192493651218106, 9.956452544367002, 11.449912743744704, 10.288421772230695]),
array([8.431678849231517, 10.901453392992671, 8.116315782793883, 8.050069058428175, 5.424809365344437, 5.49322172801956, 6.586871463912253, 4.668190623598511, 4.491513866827949, 5.111294720785227, 2.352444229841337, 2.3651836189137967, 2.333963371936645, 2.501805988839655, 1.8100031810477384, 1.3420231512314833, 1.4249384984399007, 0.37850761075592976, 0.0]),
array([10.543236661295834, 10.670663251454302]),
array([10.575290494285587, 10.047236384868398, 6.987906900679691, 5.3738848936080466]),
array([8.912838292251559, 6.21403692847767, 4.700679600659387, 4.924803007242352, 4.23123254758884, 4.638030097153396]),
array([8.867793889668471, 8.638509377601522, 7.454029927090243, 8.188655840694768, 8.036999199706717, 6.944764135462281, 5.9084485691097735, 5.6739196568842925, 5.8385141852685445, 5.049945935693616, 5.3205141570757775, 5.164060029286722, 4.887560558208756, 4.709337960597838]),
array([10.451885162096113, 10.145449683773014]),
array([10.042465900728573]),
array([9.94382593612195]),
array([7.809476068701557, 8.868236611086104, 8.8852962703692, 6.267386389907504, 4.978559874637968, 5.265347481597159, 3.8684174206812854, 3.995868051998502, 3.827435676695379, 4.548049821433009, 5.261226199713937, 3.03205870984521, 3.2218459444245595, 1.8378277975668, 2.0747280294219244, 2.4365861941489126, 2.060955455958079, 2.3067073638466558, 2.156804234844602, 2.3768848772544753, 2.212234609016627, 2.010767809296472, 1.4045127954654064, 1.1021029500296824, 0.8698952366880605, 1.734091260544155, 1.695647053130462, 1.5925941611677683, 1.2042283612949387, 1.3518357653151025, 1.7482508684128835, 0.6336729501736025, 0.6263037252226155, 0.7569656828169561, 0.27607776529789707, 0.09292986783123294, 0.0]),
array([10.197600750921909, 8.83551251688099, 9.681767860022887, 9.602671869058847, 9.969684705462114, 9.839851330846862]),
array([7.325138704138586, 9.723127983673779, 9.6362783417961]),
array([7.849962502393938, 8.81587528632126, 8.978496000838756, 9.530049947892072, 8.39731122926645, 9.304519964861116, 8.328836782111448, 9.736679663065141, 6.714653605339052]),
array([9.408405424829454, 8.601050611040241]),
array([7.741434733187303, 7.414626552631889, 8.73308484994741, 8.159169947602676, 7.2162057779453965, 5.601019433033935, 4.733102267972854, 3.9292601954403867, 3.665590963473271, 4.586060268335866, 4.544488105246264, 4.653395528264037]),
array([9.105907902777497, 9.332321813686923]),
array([7.538058147419616, 8.587617473034236, 8.256604956935291, 5.230971381564382, 4.549933361342936, 3.94791356406906, 4.1451777144974855, 4.766357770612013]),
array([6.518860331444377, 4.874262138217472, 3.9922431552757773]),
array([8.108941083156415]),
array([5.843811443868805, 5.004164590412592, 4.515985517557693, 5.003969712512245, 4.787720266116898, 2.985654720264952, 3.4689176692138353, 1.93032145075609, 1.9035023267758446, 2.4240106948092404, 2.5618799752714914, 2.4843433079341954, 2.313684296799397, 1.7487922755068173, 1.620168812555624, 1.6296884439880004]),
array([5.723232581287587, 5.600943604413235, 4.530625007939625, 4.590770455439756, 4.668393912444929, 4.202432605978718, 4.518345256755101, 2.5814101933274936, 1.9707140525058657, 0.9075951599549719, 1.5366028015875246, 0.8254291557691262, 0.9362018735735, 0.0]),
array([7.8925330211546045, 5.222097018194539, 5.2634782091145835, 2.6230748837086466, 2.7372240637593226, 2.309512004135014, 2.198671545496088, 1.684149334934156, 0.0]),
array([2.5411487170187255, 0.0]),
array([7.062717703361496, 5.715914263958712]),
array([7.879286670563875]),
array([5.9442886796236385, 5.13964146715008, 3.7790381026474504, 3.1757680576759797, 2.7457128844024634, 2.325692800874779, 1.994716328203051, 2.240854530564819, 0.9649307019886412, 0.9629656995774949, 1.4633723294673786, 0.5939766714726039, 0.0]),
array([5.938020156099044, 1.989079066714575, 2.3377770582283004, 0.0]),
array([7.190440557485331, 5.922083324047775, 6.2075388469432395, 6.826360792814582, 6.624844090640167, 6.357610184363878, 5.581963875233637, 6.379090155652856, 4.17727161162543, 4.463245281048592, 5.303311962588679, 3.986911032192565, 5.020144146256238, 5.315185920007614, 4.201646004622576, 4.689788314181972, 4.945809875556436, 4.887924867078862, 5.01088892276333, 4.221838044930637, 4.068314117834989, 3.7394922361076324, 5.060222838452524, 5.303811751065451, 3.777210796454825, 4.085138787944188, 2.8789155398325246, 1.9906777479413085, 2.4410220065454373, 1.922036409798726, 2.342457631123782, 2.517513509371404, 2.2798165445950644, 1.8713572581204518, 2.0369718401428596, 1.8312099119745773, 2.1614761242270526, 2.284215491904867, 1.8679903848031074, 1.8725667812824804, 2.4012539404044944, 2.335665385015108, 2.116188362492029, 2.333612297671764, 1.2876997949790538, 1.2293360984901918, 1.246289828279183, 1.4663327862043212, 0.9624075540677447, 1.3276528145538347, 1.2657751071404744, 0.9625986524777826, 0.8292902950591085, 0.7842172160414611, 1.701370122501345, 1.7917829855555574, 1.7480662327309457, 1.1542471323498469, 1.0635993443686758, 1.551040330349709, 1.4788253070440807, 1.3545667669994907, 1.0124359024354708, 1.5773761251175578, 1.2240721626591746, 0.9666976175186637, 1.380141740936433, 1.2084888315434315, 1.4705387231730453, 0.18055793816438404, 0.20221103172735422, 0.6415559584171149, 0.7479606440473352, 0.0]),
array([6.978566931664045]),
array([6.938880962812926, 5.997281739966459, 5.896156871595997, 3.7999881334231898, 4.313998832877531, 5.206966776526712, 3.0180855603781716, 2.833523425488748, 2.139490276174542, 2.3986153122507217, 2.433140369889221, 2.3201186671939196, 2.159264356859219, 2.4902223035389475, 1.4405402937626341, 1.0533544485569237, 1.493917992161164, 0.8898980623681219, 1.7948184539383107, 0.20983572346122203, 0.0]),
array([5.873630895159703, 6.5208202658654795, 6.983310181869642, 4.48868380915585, 3.668866833222738, 4.498560699560205, 4.420214386540167, 3.9619687906969854]),
array([6.271418546483833, 6.889941073928169, 6.817254684511638, 5.63630348821518, 4.706643838809782, 5.112463113427578, 4.393273795368196, 4.23454405988022, 4.852068191571363, 4.814159105046242, 5.121367435445928, 4.198433305327523, 3.822510413913949, 2.1443941936358746, 1.9428748550241015, 2.46417468900709, 2.366512880971374, 2.2916250989800657, 2.276633344291237, 1.1183861586810073, 1.4577538967742725, 1.571337454218511, 0.8214480075134647, 0.9759012169674063, 1.4376298882951486, 0.8358980622291163, 1.2671899988677575, 1.0145536623030642, 1.2486451015977234, 1.6347844115127104, 1.4392220217860698, 0.8888655660346766, 1.3360941304860434, 0.9121919448778905, 0.0]),
array([5.395433491706281, 6.194695729068435, 5.624393698032636, 6.128468759648447, 6.035400440826769, 4.317487324174839, 4.616629475015729, 4.136963611105148, 4.624480370630909, 5.0759369196822846, 2.868230799872907, 2.6684533399772894, 2.1430577484073647, 2.3370766144749813, 2.322944629116498, 2.277136716061461, 2.5654282978866405, 1.858002045624647, 2.230355340556964, 1.9374570421418411, 2.160600077660183, 2.0399535296690763, 2.198451001356096, 2.419165547849945]),
array([6.901573534270796, 6.6418351900377095]),
array([4.354477013984451, 5.073150899812899, 4.762507418554]),
array([6.344380010729467, 5.5193585364141375, 5.9117537516811725, 6.579943213334995, 5.562737318868761, 5.665995400135556]),
array([5.183577580633561, 5.269230963631567, 4.533361445371731, 2.7322739852411098, 1.154073971364061, 0.0]),
array([5.506717944523974, 2.7829777443529746]),
array([5.681352802525494, 3.7961916810536556, 4.638814129488589, 4.585287882140049, 3.842309860014013, 4.3396716722581745, 4.129284385363788, 4.61431943719203, 4.37576486448743, 4.513700954850657, 1.9655088220107626, 2.5428778047966722, 1.9389852147808941, 2.2553790266566502, 1.9598595995618155, 1.9796264916975317, 2.0509577999507993, 2.175292775320502]),
array([5.650016934174171, 5.230249097928067]),
array([6.128695512006625, 3.6117931233231104, 5.103699858572837, 4.90799240496122, 2.3759330425705683, 2.1324253336498553, 2.087112460103499, 0.0]),
array([2.540292322233736, 1.907513526798223, 1.4608964791706653, 1.4430712677228328]),
array([5.1189811511460945, 3.938652625603848, 5.022061055732825, 2.2176644334517417, 2.4733889423043, 2.245242687117522, 1.8660971755847293, 2.131888322737002]),
array([3.668540987164133, 4.100380187655995, 4.389623059806366, 4.840505974222532, 4.927152896363272, 2.0339299552343846, 2.253232818471663, 2.3852144562991673, 1.9867305438915077, 2.4720360466262106, 2.482574045455862, 0.9741190430971953, 0.788351674216867, 1.2161900697695587, 0.9660075397198817, 1.4713965214617428, 0.0]),
array([2.4254839647181208]),
array([3.6807456798895215]),
array([4.7834961293838365, 3.8526361450858864]),
array([3.6768978385144493, 5.225236935596438, 3.73866191460609, 3.7525233737600354, 4.158060422958452, 4.041455636357329, 2.5619220328416556, 2.259390361208009, 2.021167150512272, 2.2250427845782514, 2.3677811179430233, 1.924783565576366, 2.0510952758475556, 2.5491549627133403, 2.1045635098409714, 2.1543445557841934, 1.9530572943166131, 2.1781639093363268, 1.6236327854305581, 1.0238372675329424, 1.5803096424003433, 1.5294460237604919, 1.4035698096373053, 1.1157728071639264, 1.2367848918172881, 1.6491713340818102, 1.2596845628512314, 1.2336975927119527, 1.2326659091974155, 1.5531649363203908, 1.6673348348551147, 1.2758521910948015, 1.1826300717616538, 0.9123722348276136, 0.9408825992853418, 1.648797140227253, 0.5845059677599795, 0.17776630168090102, 0.0]),
array([2.48832589422763, 1.3707549150957656, 0.0]),
array([5.038067629188877, 3.6149879799993645, 5.0498360320595905, 4.122510829202347, 2.8757231773830294, 2.3286292816938743, 2.5572882693107966, 1.9563379742126932, 2.1661717180140565, 2.2636001768470955, 1.2154201192727319, 1.7499457558348803, 1.6604437341967553, 0.8763572118571306, 0.5591074841229361]),
array([5.171298409875613, 3.863298947797789, 4.933926833792333, 5.022428085590303, 4.210380847748518, 4.878987498642462, 5.155519249725024, 5.100855170815395, 1.9765046882016368, 2.3154456183575802, 1.8181230581349894, 1.892628147199179, 1.735430483358339, 1.5472565011146056, 0.0]),
array([4.454864897006028]),
array([4.1580009249600876, 1.8751622669011496, 1.469141939050107, 0.0]),
array([3.765291023592427, 4.190002267554081, 3.624678750427708, 1.843964618317278, 2.0384060813243807, 2.2753469924968797, 2.0897734791641542, 2.004956078507116, 2.057616068461986, 2.5445337339500003, 1.9799113769141188, 2.508101691769997, 2.0983174667311966, 2.469935028033768, 2.038924877752033, 1.6686237111236006, 1.2818744513036, 1.6667517783517045, 1.3383505480793463, 0.9412133400668337, 1.4344739346063917, 1.5273830934517523, 0.9089708017356791, 1.0390637007631223, 1.7370910408813154, 1.7446768528213321, 0.13293036687531123, 0.16279347632748298, 0.29830019941747243, 0.0]),
array([3.9709373572410214, 2.0764244917303465, 2.3211519084243433, 1.9130500639543295]),
array([3.235249513882863, 2.3012342032008877, 1.946377363830499, 2.1252621799599734, 1.5885282686091085, 0.5986370931534776, 0.0]),
array([3.95926358378067, 4.201454176178663]),
array([3.8897193748776484, 1.9410622743026082, 0.919079372657037, 1.5781721347467923, 0.7230605585031169, 0.0]),
array([3.8033382667288222]),
array([1.8040375414633416, 1.809433485215894, 2.3585479230092474, 1.8384534500281688, 1.3643384361137714, 0.8755415919814672, 0.0]),
array([3.6326794487334575, 3.735736848579615, 2.2698928260517937, 1.8531470354761144, 2.431627004836796, 0.7960628717703702, 1.72070213671286, 0.8576769797741555, 1.4985097050563698, 1.6650557187166208, 0.6998183806111871, 0.3187882162346443, 0.7094872477400849, 0.0]),
array([2.2653460483640826, 1.6962380250699007]),
array([2.2643577535138597, 0.4374290392859103, 0.0]),
array([2.2831056337245]),
array([2.9503192864303465, 2.7085144037026487, 2.1139970657259157, 2.08419881183411, 2.33256738202118, 2.196088162642196, 2.115254433361564, 2.1038470555702897, 2.416472644026087, 2.1474291938261736, 2.2368874588670775, 2.4836771448474084, 0.8465868731834103, 0.9406282430386295, 1.668054769137223, 1.6989205599881423, 1.3951324274424026, 1.2212004719539804, 1.1077260823561272, 0.2048040335322604, 0.53712128922757, 0.5297067631339225, 0.2686708664654751, 0.5428083520532748, 0.0]),
array([1.8224611683153422, 2.5637781977368093, 2.20097669970988, 2.4211749697560956, 2.0283451094374616, 1.928160101202492, 1.824243574653237, 2.2902421048689585, 1.3498702990296412, 1.0955497225396202, 1.3558180848330426, 1.741519148251279, 1.1688743193433149, 0.8308326362415658, 1.764805857410642, 0.6890139107539391, 0.28239728609319176, 0.4587915301275125, 0.4300172948617486]),
array([2.2665947928235433, 0.0]),
array([2.8837664531722194, 2.4944005051126794, 2.2882949557040995, 2.111546228555693, 0.9041403680571318, 1.3754774257197928, 1.5952927447425598, 1.2419584711473832, 1.621218235103922, 1.6370979145253273, 1.6819829734627172, 0.9815153092780656, 0.8688685182434527, 0.8975313430573962, 1.4781237193742727, 0.1980407722707942, 0.2823133921733772, 0.27073133224894186, 0.5266762474137094, 0.0]),
array([2.1379929946035365, 2.5504546442496947, 1.9444146810664713, 2.280482481927355, 2.018818533050542, 1.8613206650750462, 2.4351500486129356, 1.2905703817800733, 1.5336996360869761, 1.445321404535703, 1.4972740395543778, 0.554993651601904, 0.4028996628367145]),
array([2.388747744422331, 1.8387067316279944, 2.5327398187846923, 2.113617336602352, 1.863251109299924, 1.8058537286821061, 2.572696453820604, 0.8370991021924395, 1.3595119082737597, 0.9279549698642023, 1.5472822326769031, 1.6937075110196305, 1.713875922261618, 1.20172635476265, 0.928442929776077, 1.7023553444609372, 1.0175811877504422, 0.1948726379014205, 0.2662122708024832, 0.39752622597850573, 0.22625060463010271, 0.7074682074044291, 0.303178025423033, 0.33727219777685197, 0.18795569764536957, 0.0]),
array([1.823778664102445, 2.218249916543555, 2.4261009447068798, 1.1900195043127817, 0.9821833893840624, 1.773372636174649, 1.1906686900027355, 1.0782084671961507, 0.9192100576395541]),
array([2.055327419538813, 1.8323076816388761]),
array([1.7116252938504217, 0.0]),
array([1.9575148089504235, 2.08991852141823, 1.851079018364767, 0.8190793598186751, 1.2567936997147795, 1.312476336704837, 1.1655192185575376, 0.9766100145712281, 0.8842456096147706, 1.587041989880582, 1.1069090417265617, 0.21300607053593545, 0.25503289660375317, 0.0]),
array([2.1899251814554903, 2.1850052853989244, 1.7320044003904775, 1.2668080562974477, 1.113779052637538, 1.023199725861831, 1.7607991997853374, 0.8437021541161223, 1.5756989929480179, 1.0011460592979318, 1.4889259985827166, 1.2791649585863842, 0.41376143552949185, 0.21533045822818075, 0.5566003909288569, 0.347093411104995, 0.0]),
array([2.116224297819697, 2.108979964557299, 1.8687789055162836, 1.3067925504619, 1.4603536607511567, 1.2565834408459466, 1.6411819676557597, 0.3972372239526644, 0.48107212625926987, 0.7092040243264022, 0.6382590596079272, 0.4188395926481595, 0.0]),
array([1.7816938771916646, 0.0]),
array([1.1646497937504157, 1.5347981078605737, 0.958834482348764, 0.9550809576014149, 1.3632989134772067, 1.1599787089224, 1.3536573010678974, 0.48824786394223196, 0.31119970933784996, 0.6041521162044169, 0.0]),
array([1.805120938657981, 0.28814907348251645, 0.0]),
array([1.0703590254003523, 0.0]),
array([1.4666675465176502]),
array([1.1293177167893225, 0.9571451266538886, 0.9106093923554374, 1.355409630695547, 0.9410490762138859, 0.41815125060013597, 0.7512350258488081, 0.47403844184177635, 0.07321308679731797, 0.0]),
array([1.1654710921374982, 0.0]),
array([1.1502102885953307, 1.4378478388030038, 1.0137157042506597, 0.8411272554248084, 1.3282282242733892, 0.33149322215029675, 0.30316510780760725, 0.33416914486482846, 0.00608830498550017, 0.0]),
array([1.0076055523641587, 0.2272702357741203, 0.0]),
array([1.353567027487269, 0.524502890919238, 0.0]),
array([0.06762426391879006, 0.0]),
array([1.141456039941586, 0.8813505457078741, 1.0077762417943754, 1.1116263690163042, 0.9787155578874724, 1.1804254058697794, 1.0941986808359403, 1.0755554728984542, 0.8816881375732433, 0.9713777361261284, 1.1590015723317442, 0.1973186392902161, 0.7040036743770864, 0.2871765694795445, 0.0]),
array([1.1511438305388035, 0.9725696324453879, 0.938481512412674, 1.1632841838078294, 0.7424156423704679, 0.5298681059414829, 0.08014945405347641, 0.0]),
array([1.0012025929605943, 1.0005905716124828, 0.7997583537970764, 0.9721250796079521, 0.512147901334755, 0.0]),
array([0.6498287823931996, 0.0]),
array([0.8210430374296493, 0.6063133682199086, 0.6896404246909235]),
array([0.2019759058191321, 0.0]),
array([0.1271731740772758, 0.46198547020957137, 0.0])
]
d = [data_1]
names = ["36"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T23', 'T24', 'T25', 'T26', 'T28', 'T29', 'T30', 'T31', 'T33', 'T34', 'T35', 'T36', 'T39', 'T40', 'T41', 'T42', 'T45', 'T46', 'T47', 'T48', 'T49', 'T50', 'T51', 'T52', 'T54', 'T56', 'T57', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T66', 'T67', 'T68', 'T69', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T82', 'T83', 'T84', 'T86', 'T89', 'T91', 'T92', 'T93', 'T97', 'T99', 'T100', 'T102', 'T103', 'T104', 'T107', 'T108', 'T109', 'T111', 'T112', 'T113', 'T115', 'T117', 'T118', 'T119', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T129', 'T131', 'T133', 'T135', 'T137', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T149', 'T150', 'T151', 'T152', 'T153', 'T157', 'T158', 'T159', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T169', 'T171', 'T173', 'T175', 'T176', 'T178', 'T179', 'T180', 'T181', 'T183', 'T184', 'T186', 'T187', 'T188', 'T189', 'T192', 'T193', 'T194', 'T195', 'T197', 'T199', 'T200', 'T201', 'T203', 'T205', 'T208', 'T210']
def get_taxa_names(): return taxa_names