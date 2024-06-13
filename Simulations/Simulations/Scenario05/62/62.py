#!/usr/bin/env python
from numpy import *
data_1 = [
array([30.613928332519553, 27.032959005710737, 25.9408350245012, 27.651654206120124, 20.648113727639792, 22.5105407669761, 22.511063289354592, 20.97200950001013, 21.055206822514183, 16.766514364138043, 15.550907788977995, 15.649558369723456, 10.821900229484115, 8.791291196663893, 10.136716608804212, 3.453457215313588, 3.3611460004671962]),
array([30.85013103852556, 29.183039469296425, 27.136774002572835, 26.5065398775998, 26.398529274786725, 27.507501894650385, 24.338857050479106, 25.444818178682635, 26.580328590197528, 22.93147020474585, 21.175732693676505, 22.126418116464027, 20.444978715849622, 20.872442787631144, 22.61413762635704, 20.834835401022705, 19.10337122671337, 16.280427026081707, 16.21905756049624, 16.175567576627643, 14.443317885582434, 14.753866266191885, 12.358412351066269, 12.875146564039484, 5.095328732139041, 2.7898365660776654, 2.6706132838311203, 2.4246922217402607, 2.1500486874437144, 0.6477767652789419, 0.1480569557390723, 0.09086868927679737, 0.0]),
array([24.43135515537257, 24.782125022281, 24.67887304828293, 23.095642739492227, 23.396175742575195, 20.571044515002647, 21.20569653866806, 22.167907085311132, 16.323981533011114, 16.53625668498399, 19.441791117596605, 16.36492817956, 16.077248220667336, 14.060571995806416, 14.192427619384485, 15.34511145072295, 12.727147199148147, 12.59729979567473, 7.315572508810633, 8.330443064308882, 9.578904843792163, 11.134296033723459, 8.71257359526215, 8.400050741946764, 10.619891367506439, 7.098079812752297, 4.641838265127514, 4.25583902490306, 3.544779199274349]),
array([23.40379133365135, 24.216751270408952, 24.254483199192162, 22.81448615646485, 21.918111427487187, 22.014567836899516, 22.980518971045967, 21.25186612768288, 21.721656728729187, 20.868816023339868, 20.80271816583998, 21.19904198178075]),
array([21.047047535474412, 22.26848665375631]),
array([16.571541758435117, 16.07258495986609, 19.118806119413154, 18.352940411136682, 16.10618918161991, 17.990368455352446, 19.694810380877446, 16.450927196090213, 18.08547980181039, 16.959909837473525, 19.20075150604044, 17.2109714065685, 15.729617397509845, 14.000478167378207, 14.778898589730995, 15.103150486023619, 13.93295339371084, 14.586583995067285, 12.906655962218561, 13.808882249281018]),
array([19.281333150203952, 18.295730342940494, 17.225174535951787, 18.42529501440255, 16.72046951840584, 18.78860804903284, 18.25200736558611, 17.808188211821783, 17.88447038848209, 17.532267996301286, 15.46807886816843, 15.265090991715962, 13.95715371098549, 14.46730661410785, 14.253310538752318, 13.85529458637015, 14.377603282299411, 14.287285048007387, 14.455410384260121, 12.751206287286214, 11.561013143964503, 11.33283957288943, 11.472621260255641, 11.26803449344983, 11.037010353901163]),
array([17.867352453655858, 18.683015004083067, 16.615599383582524, 14.754545575731566, 15.160030748986772, 15.094278719414381, 15.902674622429481, 14.663750583195133]),
array([18.765120189862152, 15.413098867754414, 14.287354834243793, 15.512443417205336, 15.321403196643535, 15.606512830469454, 15.416540339936924, 14.777648801047361, 8.603088681572615, 8.279833383426613, 8.779373063917253, 8.451060664199861, 9.57912613616106, 10.260968989935629, 9.937310004107465, 11.076028479201065, 8.722937876688336, 7.9847784412945675, 7.961552255049616, 10.67211283702016, 10.027296774196484, 9.357816486897274, 6.071875608549992, 6.315965740140081, 5.34836152542459, 4.508638411314212, 3.5407127354783134, 3.1095457040007557, 3.517590999018138, 2.9171489351680324, 2.9058752172545055, 3.210566929558656, 2.8562219512128784, 1.9224092092338494, 2.340186125480248, 1.3532785861489138, 1.5016829389793052, 0.26184095573074806, 0.626313757253251, 0.5739812270527421, 0.4744892330196848, 0.0]),
array([17.66599221603094, 16.050103278782775, 15.639041488375394]),
array([16.787759546527706, 17.262287099659392, 16.758886045394707, 16.168200963132843, 15.650416738937713, 14.489664573043372, 14.516102002966571, 14.396649998220694, 14.894698471525887, 12.644194334972154, 12.26565485049915, 8.412103168481298, 7.373353432448678, 9.606093915075345, 8.788606052120151, 6.2121698463305, 4.4228896797367, 3.6086464956642845, 5.266141542415807, 3.580700916027862, 3.0508267039829797, 1.3135471056856383, 0.7648054048598181, 0.7492050941153108, 0.5398009640968482, 0.21616540108807236, 0.0]),
array([16.762001092082382, 16.545572132804796]),
array([15.451362282754681, 15.905889900793202]),
array([16.342452221854582, 14.776983805828946, 14.609864943848821, 15.834530665202305, 15.569414214236103, 13.008526548179244, 12.297550794674583, 9.920453796012923, 8.704301523648763, 10.710308556700191, 9.148179075703803, 10.424000605618978, 7.9696753628153285, 2.851406629601862, 0.723133218310223, 0.5602522894278064, 0.4455020678288231, 0.0]),
array([15.314798178833284, 14.811879613249781, 13.557975398870106, 8.185345665000565, 8.660513934481681, 10.2131473444778, 9.753880562682589, 10.994355608996454, 11.354542802168728, 9.299483363518732, 9.502925732911926, 5.45528262968212, 6.587444592202185, 5.474835461731252, 4.708737356715198, 3.417121108911925, 2.614164925736388, 3.5312620081047035, 3.4160581511840915, 3.173542747851787, 0.7782938892280404, 0.686873568359003, 0.5917490927006245, 0.6454553109314976, 0.0]),
array([16.997325411671465, 16.69975879587315, 16.348100859471934, 14.727981391703514, 14.728865578250057, 15.751637189432406]),
array([15.418586101979901, 15.176112003928898, 14.334518207187182, 13.7030210763791]),
array([15.382668596470614, 15.408884953863343, 15.496714259297805, 15.80589424981948, 15.156239319874835]),
array([15.183912523826628]),
array([14.177965591468773, 14.270801785974292, 11.983558470736297, 12.51575601261839, 13.814889020005964]),
array([12.975791305228274, 13.814550874476051, 11.233517675943983]),
array([14.99949372139543, 14.597393774726223, 15.037435472971579, 13.700693210710611]),
array([14.520070436202486]),
array([14.100119018691432, 14.489997326137773, 13.898300834510074]),
array([11.872759292293182, 12.571733134073797, 8.350052893840095, 11.430135874327767]),
array([14.0600676686739, 13.081176414267532, 12.085214079409234, 12.917621962573342, 9.660645502372548, 9.015943815765066, 8.67577689382964, 11.500627606892023, 7.48450989906216, 8.889372160217011, 7.530750890162411, 10.69341847068549, 9.22201052915517, 5.368751590328515, 5.7721345358021505, 6.049630095163956, 2.8639864817583227, 2.9995449951923203, 3.331296914543206, 2.9064951885938948, 2.328796553274475, 1.748829585722827, 0.5984500505514839, 0.725235379305821, 0.10592507978200837, 0.0]),
array([13.521540697018928, 13.530681906103467]),
array([12.453534259870954]),
array([13.048382049983209, 11.953660445355673]),
array([13.348790241175378, 13.393516283080821]),
array([12.575205566982111, 10.090340646532379, 8.9529545411737, 7.357210122470224, 8.543377193866226, 7.351364645557062, 11.46181034708289, 8.121936052688906, 9.599415864451975]),
array([13.036099649700265, 12.60674647503134]),
array([13.269821214632428, 12.764352844005893, 8.829708269515832, 10.565126037190685, 11.159838466819108, 11.30331707467017, 8.36643489652322, 9.196019025481931, 9.750559783516252, 11.100969967033686, 11.424062865517644, 11.415053649145387, 11.185501062314705, 8.51002706249513]),
array([11.832010905296697, 10.933711579810117, 11.497958731948728, 10.9806773407345, 10.678608033837865]),
array([10.266143540582435, 9.983008410038593, 10.384762633980333]),
array([12.18644894906732, 12.192871879819933, 12.73398136722985, 10.45048498849443, 9.604750855508492]),
array([12.540883491487795]),
array([12.66705709373238, 12.48301716871166]),
array([12.43630862060591, 9.805698970602196, 8.48816831357329, 9.784309312707835, 10.140766416568725, 11.4270107639015, 10.25397523765516, 9.453912681173254, 7.854453240259501, 10.684243201661149, 9.25029406013029, 9.876097503158762, 5.829225663403019, 6.9624203176746535]),
array([11.006994417459829]),
array([12.04333818590769, 9.343766989808396, 10.829789094315563, 8.932340462524259, 7.340649297472969, 10.651090255275635, 7.475057264342939, 4.750915330729333, 4.889395268997172, 3.2776707341327214, 3.0320894094242012, 3.536160768839118, 3.059426129809117, 0.675719275598592, 0.2972887374701594]),
array([11.32993638569235]),
array([11.833909726865187, 11.203709650156124, 11.234599517254544, 10.963964455700397]),
array([11.298120187355385, 11.032070879676693, 10.493203796603346]),
array([10.828709686448269, 9.096638358010464, 9.692315360431547, 8.843821755378604, 5.702756423542425, 7.1684259551799645, 5.592461457814354]),
array([10.975577854267865, 9.873251261726947, 10.677059476342833, 10.171926096205775, 10.48926865903725, 10.144632539786025]),
array([10.13050635939465, 10.540621891445493, 8.521806705163158, 5.584660865924512, 4.278394453379239, 3.5479479007621935, 1.019747034638852, 0.0]),
array([10.276229842716713]),
array([7.26248824617072, 8.873263831611137, 8.2525711267886, 8.33059337711588, 9.899362394400237, 8.194899035167198, 9.773129016888033, 6.436149720463594]),
array([9.971628059287765]),
array([11.056484651201716, 10.600907448154327, 11.063866201832372, 10.325988950270586, 8.826680906346125, 11.097285646388745, 9.226216171451538, 6.742287528770675, 7.119776367432889, 6.30188003809363, 3.600553101165019, 3.3973769521783823, 3.439400584188277, 2.949655417546718, 2.6428300011753394, 3.101320330798269, 2.9066666074664806, 2.4609576099569237, 2.064408241376049]),
array([9.228530358213918, 10.656142559991359, 7.974117425083737, 10.979287386599708, 6.4176148167983555, 5.792874485276657, 6.639881320429115, 4.164511516254986, 4.631142725784446, 3.249547200355408, 2.935835491511766, 1.6107913064416914, 0.5912919483998349, 0.7749702806486912, 0.0]),
array([9.564640557134604, 9.946572238366807, 10.086189932957915, 9.709659821019725, 10.499305190038044, 11.068096221814711, 10.978817180984013, 10.235367469834847]),
array([9.170841585517172, 8.240789748738344, 8.289037388550192, 8.660544604704826, 9.773154576538696, 9.259540257918825, 8.59061033955864, 10.423920883272084, 8.759446150389236, 9.434051563563521, 8.441398557345, 7.261560424449842, 8.379006967517164, 6.199679576561694, 5.470968274242199, 4.037535537627617, 3.401691642846908, 3.0081185887420343, 2.688922488821375, 2.0624539166236238, 2.455887469674329, 2.35738825000997, 0.43290206217124233, 0.1451507301952113, 0.3225429138258559, 0.0]),
array([7.395491121267694, 5.606543987460835, 5.07214234086925, 3.6397626066455944, 2.955325995882461, 2.5871807074695035, 2.670463740843104, 3.4767932468162286, 3.128821765090736, 2.8651326732327416, 3.42812953853445, 3.2743926226755398, 3.205574312743097, 3.206591039923215, 3.2106419001329645, 3.362169479784682, 3.1257085818350365, 3.0042728317966185, 2.710518961718663, 0.33852055422747884, 0.20005708233671937, 0.3298399254497966, 0.40255680883500605, 0.10625117373372633, 0.0]),
array([10.29430563072497, 10.40340153921107, 10.239102245571457, 10.154464112779653]),
array([9.592905078622863, 8.767417006373856, 10.01537242886248, 7.354560169180097, 7.5052428516471235, 7.019156472271649]),
array([8.786520728691158, 9.92224898607691, 9.369541199144706, 5.884140085884203, 6.7255888889820215]),
array([9.02600481927999, 5.277610606481418, 3.6400336971741765, 3.026491096405013, 3.1534149189226213, 1.9211547793638446, 0.764333865210829, 0.28883214514376054, 0.0]),
array([9.719495213950259, 9.026704973376042, 9.346286034824791, 8.722945037605832]),
array([9.039154396384324, 8.979372371632024]),
array([8.348705386054785, 7.604500469665828, 9.285794501607048, 8.550682716856713, 8.24006536083529, 7.960636981680134, 8.642494715524554]),
array([8.002696161058552, 6.338591373071204, 4.936431887407452, 4.787693006743288, 4.088762082585257, 3.5326872520336616, 2.7440267058661485, 0.6278877133019656, 0.5357097829634274, 0.5667438005150214, 0.0]),
array([8.866115381877572, 8.766235418078107, 7.7227013185469, 8.115217909003052, 8.921434693095096, 9.39629544666093, 8.43348462835804, 9.41003913077111, 7.202931317244924, 5.8539384720592, 5.602543422007432, 5.089700866382696, 4.889061219180581, 4.681265685584224, 2.740269571360038, 3.5124837508883067, 2.9706915833655567, 3.4005414153212574, 2.8552758106447964, 3.0715570529425977, 3.3555730169947666, 3.196687733620679, 1.8120065573562263, 2.543578773224632, 1.913040091476374, 1.7219242921663669]),
array([9.065674556705469]),
array([9.113187405159414, 7.625067924452466, 7.911783997748006, 8.85670844078798, 8.289480264276564, 9.011019282030551, 8.320665950398888, 9.162056406788903, 8.849966495959814, 7.673755492308562, 7.74050440465242, 7.680880901596069]),
array([8.424187589976187, 8.215021413564493, 6.8694923400985655]),
array([8.833147621665702]),
array([9.007791622432888, 8.935166841546671]),
array([8.408018433506244, 7.976106867940509, 8.428997848975555, 7.825914236447675, 8.624482229664633, 7.145550330987815]),
array([8.949077960451978, 7.4613944584676855, 8.37695694501047, 8.694817625575263, 3.339289786594597, 3.0524018377793816, 3.066503137700253]),
array([8.043607712513394, 8.627167333666911, 8.85810092819407, 6.230247821845176, 6.965922187428688, 6.853299712866415, 6.3691139383841975, 6.789153539960269, 4.983069292112048, 4.2700986190602, 4.505497851211175, 3.9070120398128854, 4.941824236753255, 2.918732541908792, 3.1023246467015446, 3.3886732809004934, 3.253525264394539, 3.3423128831100817, 3.3569924664961865, 3.0216821259673137, 2.735654383843212, 3.0589621571763885, 2.876903662363012, 2.6810223240007387, 2.3872172146862294, 1.9097905695953736, 1.9719463422906407, 1.5594144830023455, 0.1873924698636279, 0.3404044338089886, 0.7174574112445035, 0.24571631449582243, 0.5098233710243498, 0.450866771971188, 0.7059619177582029, 0.5895092313593802, 0.32192621705950064, 0.07732897274192352, 0.0]),
array([8.187695040045366]),
array([7.4225836636253835, 5.990661816430122, 5.631299529694974, 6.590216962383446, 4.456935736047386, 3.523884743587069, 3.285692477690806, 2.9058944178029003, 3.3250838991672347, 3.0674563489197504, 3.5336879381025494, 2.6901375916791777, 3.5344267321196243, 2.801662884523574, 2.9607864474428722, 3.4450783857788627]),
array([7.544227199216978, 7.292154097154064, 8.459406907910276, 8.2018377960804, 7.942964519462445, 7.367682350835386, 7.656262537148049, 5.351537585470544, 6.015712363683973, 7.171866126636916, 6.601683642361236, 4.080288405404694, 3.831057382251603, 4.739027557617821, 3.3432693780740834, 3.5318310693271755, 3.066746900868825, 2.2115351762223865]),
array([7.84229364151465, 8.455549616402447, 8.269161814009362]),
array([7.759693456535566, 7.793332108542873, 7.031358682419967]),
array([7.9509756562052605, 7.926345917041779, 7.402427817356134, 7.799925821002705, 7.324363890185245, 7.9870980121818675, 6.918426934936544, 7.078230506784718]),
array([7.960940581554711, 6.699494000507464, 6.396610225284989, 3.553048954907217, 3.337472899246926, 2.831537787226486, 2.610030116516934, 3.3657011317505656, 3.1328337008739364, 2.4476015914955274, 0.7610894532369882, 0.7709279651858366, 0.7767312580036646, 0.08451936535806917, 0.0]),
array([6.825420436011773]),
array([7.527867835247992, 7.717043515339923, 7.486107694637425]),
array([7.556474670819465, 7.658338856828083]),
array([7.326182852454742, 6.959732788526382]),
array([7.333802362872316, 7.129357459594544, 3.8685327011947814, 5.08010990581225, 2.7023758102674282, 2.851810426566364, 3.4665203613154807, 0.7752020096675236, 0.6103766247585842, 0.452259054086098, 0.40468561370998113, 0.655047008217237, 0.17286341496216695, 0.0]),
array([6.977156113350036, 5.501436909551356, 4.859780571220961]),
array([7.465076216874963, 7.519705099539986, 7.348489613070264, 5.3464819124460945, 6.87854470249761, 5.098968100303867, 4.680707959864427, 4.728790292720919, 4.867499875633239]),
array([6.04052203210056, 5.4933240870348765, 6.509662420675751, 3.2955745434310133, 3.4618589505296953, 2.9446632945692706, 2.991068704389369, 3.586187236227849, 3.5209655539366866, 2.7812463318295193, 3.0339329402968787, 2.60095841260812, 1.9828324772505852, 1.805959041111335, 2.265414806504788]),
array([7.147858597438291]),
array([7.247113294881243]),
array([6.81385560703112]),
array([5.258986488660524, 5.111821988635374, 2.859888959766885, 2.8760760919878723, 2.4881077220101466, 0.3840168121980725, 0.4131595301053355, 0.04687404428339942, 0.0]),
array([4.5737087432490915, 2.959657920166484, 2.4484100950660626, 1.264078326975043, 0.23549554376436677, 0.5650498024046069, 0.0]),
array([4.4745907770181175, 3.00091232679095, 3.314860047983829, 3.1161527591270124, 2.3463353500627853, 2.529465227373398, 2.495655814119582, 0.0]),
array([6.28332257234139, 5.072204693113051]),
array([6.629405037185309, 0.0]),
array([6.084473956156758, 5.773696858255896]),
array([5.5041368674835685, 6.4446267354234275, 5.637264678303829, 6.317694010658942, 6.040296928343888, 4.1683979263152775, 3.5467812632729583, 3.440420568457222, 3.5937837887974986, 3.4023594728392963, 3.142767300253685, 3.2548523259610183]),
array([6.104099343655371, 4.099818291346183, 3.2087172184437125, 3.376594611244565, 3.4441129732705917, 2.4003029852540387, 2.5134619121405644, 2.561659337776882]),
array([5.7436207267596835, 6.421176421339118]),
array([4.827774982315047, 3.591722774081056]),
array([5.708175927428723, 4.620614035986407, 2.589746373659366, 1.937347831938756, 1.9023395827813336, 0.32570961730841375, 0.43901674898968374, 0.6247016817047183, 0.0]),
array([5.739605216518494, 3.4951187161973967, 0.7333259287982898, 0.4915495594675744, 0.0]),
array([5.268136538467557]),
array([5.528133189096485, 4.426746261151956, 4.683049307337107]),
array([5.48807160485841]),
array([4.075271521760849, 4.783279750053684, 3.438337494901215, 2.861721796516955, 3.107735295701927, 3.120624938757046, 3.483620972750201, 2.6930791200080795, 3.54826988092053, 3.0959264393032497, 1.9308579043317007, 2.3232745549125156, 2.343144421685836]),
array([5.224572592866789, 2.8513651211348554, 3.5950730651686094, 3.424498743546543, 3.338677295499802, 2.48478400698195, 2.187587802357072, 1.9184515331013459, 0.58654441029328, 0.7349755512071602, 0.0]),
array([2.6710480102227896, 3.550918390282326]),
array([2.992992427000591, 2.7942520687054544, 2.6790140847195207, 2.6700619164081987, 2.8431588426055057, 3.493185936158915, 3.385223578536979, 2.236416707539055, 2.5378741998155863, 0.721614847820632, 0.2904685589657202, 0.5626050764745878, 0.3845216968985496]),
array([2.6727153555338883, 2.878271646113194]),
array([3.778623663191168, 4.297974482860861, 3.0066736234934512, 3.164917919254063, 3.2996218991721733, 2.366800568238671, 2.5533159867831077]),
array([4.685238581976519, 3.678404759532314, 3.3408371562467982, 3.085668986918201, 3.233373869761292, 3.366728222427039, 2.986906545749396, 2.912063266967211, 0.20027568058628453, 0.1344736845212825, 0.3939295249489775, 0.31443127273216404, 0.021748253168378723, 0.002881019294592696, 0.0]),
array([2.808237914957319, 3.376344884286696, 2.6442098058140204, 3.1879297462475664, 2.0401741504934647, 2.563085595233667, 0.8922116963755741, 0.6311495559255331, 0.2087288860606069, 0.4235684798297281, 0.2585345500900589, 0.0]),
array([4.306110259964493, 4.3342921789154865, 2.8509187643156246, 2.9971351005155675, 2.7452696554720877, 2.862058742269449, 2.354445010714359]),
array([2.9435215585414793]),
array([3.929679566446529, 2.9751074417912733, 3.3298625297855855, 3.5569993176497654, 3.1841877435566857, 2.341688388296233, 1.8688890874958348, 0.10651892116113326, 0.0]),
array([3.857037475226993, 4.378119530141845, 3.9315324357560586, 3.355163290333324, 3.2544929660383715, 0.6589968453013683, 0.6458846272723197]),
array([3.8831182551221715, 3.7232737421429123, 4.310307583065586, 2.9641157986908877, 3.550939083924283, 3.151533838186766, 3.4200484372677527, 2.839482908667093, 3.0653534177944763, 2.5711191739892274, 2.2990133154113717, 0.6753435070987386, 0.4580164894247381, 0.5500736866819911, 0.07393607996633125, 0.0]),
array([3.1828366185966077, 2.738433675643459, 3.4924599010206547, 3.168535772119263, 3.0150021474146333, 2.680591649163365, 3.010530255146561, 3.597271325448113, 3.1612622009609574, 2.726583958702391, 2.7035618662377505, 2.9360525042836274, 2.357590903619151, 2.5620502286337743, 1.529468293116963, 0.583131090940316, 0.4037567558925306, 0.2691301454520314, 0.37415241304675584, 0.5670293889902122, 0.5841159636004618, 0.20936663480794293, 0.6101235098903328, 0.12985435363972975, 0.5125896075024585, 0.044238146075712415, 0.0]),
array([4.244486441517564]),
array([3.392098439633671, 3.1461976254933144, 0.3355659599838158, 0.612104244449184, 0.0]),
array([3.886830858090654, 2.6196635563744546, 2.6442859438719846, 2.626902584905393, 3.474707832419378, 3.121967525699821, 3.340530846861527, 3.365718940577758]),
array([3.1994355988357293]),
array([2.5976349165218147, 2.918744319788642, 2.8194510115367057, 2.855190638941903, 2.453638493613557, 0.41857805575959617, 0.3275937645457364, 0.3477417674334812, 0.2996260305763066, 0.10418801817530188, 0.0]),
array([3.0930143621050457, 3.2015770802407495]),
array([3.256292261942657, 2.86284024310183, 2.8599972712660913, 3.344727668560178, 2.413220279738509, 1.835089993671157, 2.5739657871401134, 0.7595811125501433, 0.335868274155578, 0.3110230972810962, 0.2827611690992988, 0.09052347094731374, 0.0]),
array([3.310645782627126, 3.241768325980288, 3.1668034407147108, 3.271497893701083, 2.6248701691704395, 2.4438741246089437, 2.452492144364694, 0.9890213740423006, 1.403218402772308, 0.5599298054493433, 0.3604589359944445, 0.0]),
array([3.243040217324987, 2.7809287673771586, 3.5515128025514513, 2.9065367537689353, 3.27611982473524, 2.5794619845940923, 0.5477931118288956, 0.25554831741533, 0.7415588718294747, 0.0]),
array([2.964808508696182, 3.510171341905374, 2.8351219262074494, 3.5419894867343342, 2.654885993857203, 3.2594421330870964, 3.2651678792629895, 3.5801996771718705, 3.000965114649152, 3.269479748621519, 3.5491468154050714, 2.0089445093153664, 0.17627750194311076, 0.20061175787721386, 0.0]),
array([3.4893308601301256]),
array([2.040340711394969, 1.981670904750325, 0.0]),
array([3.197333369522198, 2.8775874182744388, 3.0449315469233063, 2.7273106777272673, 3.055632234145013, 2.916472627433505, 2.9849691268799297, 2.743874656775828]),
array([3.0676134811202487, 3.1030441455699327, 2.802966821924109, 3.1263025710122974, 3.137351072385765, 2.8619736379428207, 2.704582121697798, 3.2316646909278917, 2.946131530089038, 1.953580374378092, 1.576153413946491, 0.48957001634802666, 0.7012492324201411, 0.2980510062996792, 0.2620873031488198, 0.4790107851804981, 0.0]),
array([2.2564291364786033, 2.4704444622018835, 2.4770938667893603, 1.87174503706312, 1.8570480542230565]),
array([2.771263262629017, 2.8492878132755806, 2.7951885895458135, 2.329141636794782, 0.9511133992493548, 0.39881097848001223, 0.450551891515341, 0.4417411574816633, 0.0]),
array([2.6901969018656953, 2.88252340415684, 2.8505250840279346]),
array([2.6145971248010413, 2.494843349145571, 1.9064855994549967, 2.01477822100276]),
array([2.7559498470513244]),
array([2.9144923695420486, 2.6483548925572036, 2.8177332742828405]),
array([2.0093377844256977, 0.2709372887386594, 0.40523203082815806, 0.0]),
array([2.7127523538199765, 2.7872822546949156, 2.1967058102523653]),
array([2.2023122201404566, 0.14723263595777247, 0.4415847646242978, 0.6117625733191118, 0.7096448325512811, 0.7566564927820533, 0.16368561810654925, 0.0]),
array([0.8497085166389343, 0.6184568050321242]),
array([1.9409665053521437, 2.403885850201325, 0.3283686018540698]),
array([2.506809046000723, 2.2920269892375056, 1.2214166507684667, 0.20147743884658575, 0.6238713981947196, 0.16139784663845058, 0.20671699339148808, 0.7205105474210859, 0.49008558052912693]),
array([1.9916284334275922, 1.2136290850853682, 0.639811007102155]),
array([2.1006230291050567, 2.054807323935294, 0.9255441853795424, 0.7153324133565345, 0.20936670016081171, 0.6367847355452289, 0.6430350002278655, 0.18591123100138673, 0.0]),
array([2.200287687818172, 2.295659985251341, 0.9898645991346234, 0.5480383220369835, 0.729801968223682, 0.24116829985638077, 0.3749411859333721, 0.0]),
array([2.179161263514195, 2.1791738295300482, 0.39346703012251627, 0.3229030706019143, 0.5627121069887244, 0.2034142640206561, 0.05931703189088566, 0.0]),
array([2.0625431304928337, 1.306841667494867, 0.67240144920712, 0.7505079875724264, 0.7140769029502886, 0.38259559663300635, 0.03208877465872248, 0.0]),
array([1.8498351556370234, 2.1119760138771664, 2.1130581288345094]),
array([1.873362181889557, 1.8754766779201806, 1.8140146881602237, 1.7684576760062591, 0.5066992304802891, 0.6851346394175207, 0.7412839300470534, 0.7418913916250524, 0.3091800148890945, 0.5944968034958975, 0.0]),
array([1.9500167426969925, 1.8170395872055747]),
array([1.8559269258731994, 1.5102855619509723, 0.32612104067450165, 0.0]),
array([1.7845469084123948]),
array([0.3939169023881548, 0.6222354852202427, 0.1726059782528846, 0.592232276282409, 0.0]),
array([0.6924279631166197, 0.40683507083448495, 0.40387763796684456, 0.2783035460801737, 0.2979285643925377, 0.006107089059844875, 0.0]),
array([0.7244774125394904, 0.6905065779734625]),
array([1.083720733671607, 0.581840123586711, 0.5761804129745357, 0.5743918575724867, 0.7444064665389138, 0.6355578886699212, 0.4489041559368333, 0.0]),
array([0.20680496860387054, 0.12678156138392616, 0.0]),
array([0.6234165696337114]),
array([0.45822655670004026, 0.0]),
array([0.24390045414148342, 0.45387034168673235, 0.34324549334625015, 0.0]),
array([0.7464885446416782, 0.4501718085710272, 0.0]),
array([0.17972383457987418, 0.5551059083642906, 0.10355145749864082, 0.0]),
array([0.40201925871483984, 0.7637099016969818, 0.6695285574536903, 0.3607080948986526, 0.0]),
array([0.3842359679278958, 0.0]),
array([0.5294456670289693]),
array([0.4212393326973084, 0.6255899114307597, 0.6427157767390947, 0.32218524344040145, 0.40313906506240865, 0.04454797841550043, 0.0]),
array([0.29050250854353743, 0.3855946596612571, 0.4360957307155252, 0.0]),
array([0.4163068615118231]),
array([0.467763432681376, 0.47415456414206025, 0.2171188053674415, 0.427089441678629, 0.3717479280729161, 0.048295414096023365, 0.0]),
array([0.4335268258232179, 0.1936755009686345, 0.0]),
array([0.3476189575018067, 0.3369594431416755, 0.0]),
array([0.013241286289440015, 0.0])
]
d = [data_1]
names = ["62"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T28', 'T29', 'T30', 'T31', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T41', 'T42', 'T46', 'T47', 'T48', 'T49', 'T51', 'T52', 'T53', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T66', 'T68', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T86', 'T87', 'T90', 'T92', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T109', 'T111', 'T114', 'T115', 'T116', 'T117', 'T118', 'T119', 'T121', 'T124', 'T125', 'T126', 'T127', 'T129', 'T130', 'T131', 'T133', 'T134', 'T135', 'T137', 'T138', 'T139', 'T140', 'T142', 'T144', 'T145', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T168', 'T169', 'T170', 'T171', 'T172', 'T173', 'T175', 'T176', 'T177', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'T186', 'T188', 'T189', 'T190', 'T192', 'T199', 'T201', 'T204', 'T205', 'T208', 'T209', 'T210', 'T211', 'T212', 'T215', 'T216', 'T217', 'T218', 'T219', 'T221', 'T222', 'T226']
def get_taxa_names(): return taxa_names