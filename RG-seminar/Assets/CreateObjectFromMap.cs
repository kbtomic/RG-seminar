using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;

public class CreateObjectFromMap : MonoBehaviour
{
    //43.51293893730756 (ULat), 16.434653728432924 (ULong) gornja liva
    //43.51175528819933(LLAT), 16.435623998944898(LLong) donja desna
    public class Map
    {
        public float Width { get; set; }
        public float Height { get; set; }
        public float UpperEdgeDiagonalLatitude { get; set; }
        public float UpperEdgeDiagonalLongitude { get; set; }
        public float LowerEdgeDiagonalLatitude { get; set; }
        public float LowerEdgeDiagonalLongitude { get; set; }

        public void CalculateMapWidthAndHeight(float upperEdgeLatitude, float upperEdgeLongitude, float lowerEdgeLatitude, float lowerEdgeLongitude)
        {
            this.UpperEdgeDiagonalLatitude = upperEdgeLatitude;
            this.UpperEdgeDiagonalLongitude = upperEdgeLongitude;
            this.LowerEdgeDiagonalLatitude = lowerEdgeLatitude;
            this.LowerEdgeDiagonalLongitude = lowerEdgeLongitude;

            GeoUTMConverter UpperEdgeDiagonal = new GeoUTMConverter();
            GeoUTMConverter LowerEdgeDiagonal = new GeoUTMConverter();

            UpperEdgeDiagonal.ToUTM(upperEdgeLatitude, upperEdgeLongitude);
            LowerEdgeDiagonal.ToUTM(lowerEdgeLatitude, lowerEdgeLongitude);

            if (LowerEdgeDiagonal.Zone == UpperEdgeDiagonal.Zone)
                this.Width = (float)Math.Abs(UpperEdgeDiagonal.X - LowerEdgeDiagonal.X);
            else
                this.Width = (float)Math.Abs((1000000 - UpperEdgeDiagonal.X) - LowerEdgeDiagonal.X);

            if (upperEdgeLatitude > 0 && lowerEdgeLatitude < 0) //if map is both on north and south hemi
                this.Height = (float)Math.Abs(UpperEdgeDiagonal.Y - (10000000 - LowerEdgeDiagonal.Y));
            else
                this.Height = (float)Math.Abs(UpperEdgeDiagonal.Y - LowerEdgeDiagonal.Y);
        }
    }

    public class ObjectMap
    {
        public float LeftLowerEdgeX { get; set; }
        public float LeftLowerEdgeY { get; set; }
        public float TranslationOfLeftLowerEdgeX { get; set; }
        public float TranslationOfLeftLowerEdgeY { get; set; }

        public void CreateObjectMap(Map map)
        {
            GameObject objectMap = GameObject.CreatePrimitive(PrimitiveType.Cube);
            objectMap.name = "map";

            objectMap.transform.localScale += new Vector3(map.Width - 1, 0, map.Height - 1); //we want object to be in xz level, so real y coordinate iz z in unity

            float initialLowerEdgeX = objectMap.transform.position.x;
            float initialLowerEdgeY = objectMap.transform.position.z;

            this.TranslateMapInTheOriginOfCoordinateSystem(objectMap, map);

            this.CalculateTranslationXY(objectMap, initialLowerEdgeX, initialLowerEdgeY);

            //put real map as texture on game object map
            Texture2D texture = new Texture2D(128, 128);
            var fileData = File.ReadAllBytes("C:/StariPlac.png");
            texture.LoadImage(fileData);
            Color[] pix = texture.GetPixels(); // get pixel colors
            for (int i = 0; i < pix.Length; i++)
                pix[i].a = pix[i].grayscale; // set the alpha of each pixel to the grayscale value
            texture.SetPixels(pix); // set changed pixel alphas
            texture.Apply(); // upload texture to GPU
            objectMap.GetComponent<Renderer>().material.mainTexture = texture;


        }

        //public void PutMapAsMaterial(GameObject objectMap)
        //{
        //    Texture2D tex;
        //    tex = new Texture2D(4, 4, TextureFormat.DXT1, false);
        //    string url = "Pic"
        //    using (WWW www = new WWW(url))
        //    {
        //        yield return www;
        //        www.LoadImageIntoTexture(tex);
        //        objectMap.GetComponent<Renderer>().material.mainTexture = tex;
        //    }
        //}

        //public IEnumerator LoadTextureFromCache(string filePath, GameObject objectMap) //ne radi
        //{
        //    if (!File.Exists(filePath))
        //    {
        //        yield break;
        //    }
        //    var www = UnityWebRequestTexture.GetTexture("file://" + filePath);
        //    yield return www.SendWebRequest();
        //    //texture loaded
        //    var texture = DownloadHandlerTexture.GetContent(www);
        //    objectMap.GetComponent<Renderer>().material.mainTexture = texture;

        //}

        public void TranslateMapInTheOriginOfCoordinateSystem(GameObject objectMap, Map map)
        {
            objectMap.transform.Translate(map.Width / 2, 0, map.Height / 2);
            this.LeftLowerEdgeX = objectMap.transform.position.x;
            this.LeftLowerEdgeY = objectMap.transform.position.z;
        }

        public void CalculateTranslationXY(GameObject objectMap, float initialLowerEdgeX, float initialLowerEdgeY)
        {
            this.TranslationOfLeftLowerEdgeX = Math.Abs(initialLowerEdgeX - objectMap.transform.position.x);
            this.TranslationOfLeftLowerEdgeY = Math.Abs(initialLowerEdgeY - objectMap.transform.position.z);
        }
    }

    public class GeoUTMConverter
    {
        public double Latitude { get; set; }
        public double Longitude { get; set; }
        public double Zone { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public Hemisphere Hemi { get; set; }

        private double pi = 3.14159265358979;
        private double sm_a = 6378137.0;
        private double sm_b = 6356752.314;
        //    private double sm_EccSquared = 6.69437999013e-03;
        private double UTMScaleFactor = 0.9996;

        public enum Hemisphere
        {
            Northern = 0,
            Southern = 1
        }

        public GeoUTMConverter()
        {

        }

        public void ToUTM(double Latitude, double Longitude)
        {
            this.Latitude = Latitude;
            this.Longitude = Longitude;

            Zone = Math.Floor((Longitude + 180.0) / 6) + 1;
            GeoUTMConverterXY(DegToRad(Latitude), DegToRad(Longitude), Zone);
        }

        public void ToLatLon(double X, double Y, int zone, Hemisphere Hemi)
        {
            this.X = X;
            this.Y = Y;

            this.Zone = zone;

            if (Hemi == Hemisphere.Northern)
            {
                UTMXYToLatLon(X, Y, false);
            }
            else
            {
                UTMXYToLatLon(X, Y, true);
            }
        }

        private double DegToRad(double degrees)
        {
            return (degrees / 180.0 * pi);
        }

        private double RadToDeg(double radians)
        {
            return (radians / pi * 180.0);
        }

        private double MetersToFeet(double meters)
        {
            return (meters * 3.28084);
        }

        private double FeetToMeters(double feet)
        {
            return (feet / 3.28084);
        }

        private double ArcLengthOfMeridian(double phi)
        {
            double alpha, beta, gamma, delta, epsilon, n;
            double result;

            /* Precalculate n */
            n = (sm_a - sm_b) / (sm_a + sm_b);

            /* Precalculate alpha */
            alpha = ((sm_a + sm_b) / 2.0)
               * (1.0 + (Math.Pow(n, 2.0) / 4.0) + (Math.Pow(n, 4.0) / 64.0));

            /* Precalculate beta */
            beta = (-3.0 * n / 2.0) + (9.0 * Math.Pow(n, 3.0) / 16.0)
               + (-3.0 * Math.Pow(n, 5.0) / 32.0);

            /* Precalculate gamma */
            gamma = (15.0 * Math.Pow(n, 2.0) / 16.0)
                + (-15.0 * Math.Pow(n, 4.0) / 32.0);

            /* Precalculate delta */
            delta = (-35.0 * Math.Pow(n, 3.0) / 48.0)
                + (105.0 * Math.Pow(n, 5.0) / 256.0);

            /* Precalculate epsilon */
            epsilon = (315.0 * Math.Pow(n, 4.0) / 512.0);

            /* Now calculate the sum of the series and return */
            result = alpha
                * (phi + (beta * Math.Sin(2.0 * phi))
                    + (gamma * Math.Sin(4.0 * phi))
                    + (delta * Math.Sin(6.0 * phi))
                    + (epsilon * Math.Sin(8.0 * phi)));

            return result;

        }

        private double UTMCentralMeridian(double zone)
        {
            return (DegToRad(-183.0 + (zone * 6.0)));
        }

        private double FootpointLatitude(double y)
        {
            double y_, alpha_, beta_, gamma_, delta_, epsilon_, n;
            double result;

            /* Precalculate n (Eq. 10.18) */
            n = (sm_a - sm_b) / (sm_a + sm_b);

            /* Precalculate alpha_ (Eq. 10.22) */
            /* (Same as alpha in Eq. 10.17) */
            alpha_ = ((sm_a + sm_b) / 2.0)
                * (1 + (Math.Pow(n, 2.0) / 4) + (Math.Pow(n, 4.0) / 64));

            /* Precalculate y_ (Eq. 10.23) */
            y_ = y / alpha_;

            /* Precalculate beta_ (Eq. 10.22) */
            beta_ = (3.0 * n / 2.0) + (-27.0 * Math.Pow(n, 3.0) / 32.0)
                + (269.0 * Math.Pow(n, 5.0) / 512.0);

            /* Precalculate gamma_ (Eq. 10.22) */
            gamma_ = (21.0 * Math.Pow(n, 2.0) / 16.0)
                + (-55.0 * Math.Pow(n, 4.0) / 32.0);

            /* Precalculate delta_ (Eq. 10.22) */
            delta_ = (151.0 * Math.Pow(n, 3.0) / 96.0)
                + (-417.0 * Math.Pow(n, 5.0) / 128.0);

            /* Precalculate epsilon_ (Eq. 10.22) */
            epsilon_ = (1097.0 * Math.Pow(n, 4.0) / 512.0);

            /* Now calculate the sum of the series (Eq. 10.21) */
            result = y_ + (beta_ * Math.Sin(2.0 * y_))
                + (gamma_ * Math.Sin(4.0 * y_))
                + (delta_ * Math.Sin(6.0 * y_))
                + (epsilon_ * Math.Sin(8.0 * y_));

            return result;
        }

        private double[] MapLatLonToXY(double phi, double lambda, double lambda0)
        {
            double[] xy = new double[2];

            double N, nu2, ep2, t, t2, l;

            double l3coef, l4coef, l5coef, l6coef, l7coef, l8coef;

            //      double tmp;



            /* Precalculate ep2 */

            ep2 = (Math.Pow(sm_a, 2.0) - Math.Pow(sm_b, 2.0)) / Math.Pow(sm_b, 2.0);



            /* Precalculate nu2 */

            nu2 = ep2 * Math.Pow(Math.Cos(phi), 2.0);



            /* Precalculate N */

            N = Math.Pow(sm_a, 2.0) / (sm_b * Math.Sqrt(1 + nu2));



            /* Precalculate t */

            t = Math.Tan(phi);

            t2 = t * t;
            //        tmp = (t2 * t2 * t2) - Math.Pow(t, 6.0);

            /* Precalculate l */
            l = lambda - lambda0;

            /* Precalculate coefficients for l**n in the equations below
               so a normal human being can read the expressions for easting
               and northing
               -- l**1 and l**2 have coefficients of 1.0 */

            l3coef = 1.0 - t2 + nu2;

            l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);

            l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2
                - 58.0 * t2 * nu2;

            l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2
                - 330.0 * t2 * nu2;

            l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
            l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);

            /* Calculate easting (x) */
            xy[0] = N * Math.Cos(phi) * l
                + (N / 6.0 * Math.Pow(Math.Cos(phi), 3.0) * l3coef * Math.Pow(l, 3.0))
                + (N / 120.0 * Math.Pow(Math.Cos(phi), 5.0) * l5coef * Math.Pow(l, 5.0))
                + (N / 5040.0 * Math.Pow(Math.Cos(phi), 7.0) * l7coef * Math.Pow(l, 7.0));

            /* Calculate northing (y) */
            xy[1] = ArcLengthOfMeridian(phi)
                + (t / 2.0 * N * Math.Pow(Math.Cos(phi), 2.0) * Math.Pow(l, 2.0))
                + (t / 24.0 * N * Math.Pow(Math.Cos(phi), 4.0) * l4coef * Math.Pow(l, 4.0))
                + (t / 720.0 * N * Math.Pow(Math.Cos(phi), 6.0) * l6coef * Math.Pow(l, 6.0))
                + (t / 40320.0 * N * Math.Pow(Math.Cos(phi), 8.0) * l8coef * Math.Pow(l, 8.0));


            return xy;
        }

        private double[] MapXYToLatLon(double x, double y, double lambda0)
        {
            double[] latlon = new double[2];

            double phif, Nf, Nfpow, nuf2, ep2, tf, tf2, tf4, cf;
            double x1frac, x2frac, x3frac, x4frac, x5frac, x6frac, x7frac, x8frac;
            double x2poly, x3poly, x4poly, x5poly, x6poly, x7poly, x8poly;

            /* Get the value of phif, the footpoint latitude. */
            phif = FootpointLatitude(y);

            /* Precalculate ep2 */
            ep2 = (Math.Pow(sm_a, 2.0) - Math.Pow(sm_b, 2.0))
                  / Math.Pow(sm_b, 2.0);

            /* Precalculate cos (phif) */
            cf = Math.Cos(phif);

            /* Precalculate nuf2 */
            nuf2 = ep2 * Math.Pow(cf, 2.0);

            /* Precalculate Nf and initialize Nfpow */
            Nf = Math.Pow(sm_a, 2.0) / (sm_b * Math.Sqrt(1 + nuf2));
            Nfpow = Nf;

            /* Precalculate tf */
            tf = Math.Tan(phif);
            tf2 = tf * tf;
            tf4 = tf2 * tf2;

            /* Precalculate fractional coefficients for x**n in the equations
               below to simplify the expressions for latitude and longitude. */
            x1frac = 1.0 / (Nfpow * cf);

            Nfpow *= Nf;   /* now equals Nf**2) */
            x2frac = tf / (2.0 * Nfpow);

            Nfpow *= Nf;   /* now equals Nf**3) */
            x3frac = 1.0 / (6.0 * Nfpow * cf);

            Nfpow *= Nf;   /* now equals Nf**4) */
            x4frac = tf / (24.0 * Nfpow);

            Nfpow *= Nf;   /* now equals Nf**5) */
            x5frac = 1.0 / (120.0 * Nfpow * cf);

            Nfpow *= Nf;   /* now equals Nf**6) */
            x6frac = tf / (720.0 * Nfpow);

            Nfpow *= Nf;   /* now equals Nf**7) */
            x7frac = 1.0 / (5040.0 * Nfpow * cf);

            Nfpow *= Nf;   /* now equals Nf**8) */
            x8frac = tf / (40320.0 * Nfpow);

            /* Precalculate polynomial coefficients for x**n.
               -- x**1 does not have a polynomial coefficient. */
            x2poly = -1.0 - nuf2;

            x3poly = -1.0 - 2 * tf2 - nuf2;

            x4poly = 5.0 + 3.0 * tf2 + 6.0 * nuf2 - 6.0 * tf2 * nuf2
                - 3.0 * (nuf2 * nuf2) - 9.0 * tf2 * (nuf2 * nuf2);

            x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2;

            x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2
                + 162.0 * tf2 * nuf2;

            x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2);

            x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2);

            /* Calculate latitude */
            latlon[0] = phif + x2frac * x2poly * (x * x)
                + x4frac * x4poly * Math.Pow(x, 4.0)
                + x6frac * x6poly * Math.Pow(x, 6.0)
                + x8frac * x8poly * Math.Pow(x, 8.0);

            /* Calculate longitude */
            latlon[1] = lambda0 + x1frac * x
                + x3frac * x3poly * Math.Pow(x, 3.0)
                + x5frac * x5poly * Math.Pow(x, 5.0)
                + x7frac * x7poly * Math.Pow(x, 7.0);

            return latlon;
        }

        private void GeoUTMConverterXY(double lat, double lon, double zone)
        {
            double[] xy = MapLatLonToXY(lat, lon, UTMCentralMeridian(zone));

            xy[0] = xy[0] * UTMScaleFactor + 500000.0;
            xy[1] = xy[1] * UTMScaleFactor;
            if (xy[1] < 0.0)
                xy[1] = xy[1] + 10000000.0;

            this.X = xy[0];
            this.Y = xy[1];

            //this.X = FeetToMeters(this.X);
            //this.Y = FeetToMeters(this.Y);
        }

        private void UTMXYToLatLon(double x, double y, bool southhemi)
        {
            double cmeridian;

            x -= 500000.0;
            x /= UTMScaleFactor;

            /* If in southern hemisphere, adjust y accordingly. */
            if (southhemi)
                y -= 10000000.0;

            y /= UTMScaleFactor;

            cmeridian = UTMCentralMeridian(Zone);
            double[] latlon = MapXYToLatLon(x, y, cmeridian);

            this.Latitude = RadToDeg(latlon[0]);
            this.Longitude = RadToDeg(latlon[1]);


        }

    }
    // Start is called before the first frame update
    void Start()
    {
        //43.51293893730756 (ULat), 16.434653728432924 (ULong) gornja liva
        //43.51175528819933(LLAT), 16.435623998944898(LLong) donja desna

        //0.0009978417767449134, 9.374125157568141 sjeverna i juzna polutka
        //-0.0007369909727086933, 9.377088232797714

        //-0.006765411712556136, 9.367515224307253 ove su obe juzne polutke 
        //-0.009284593284594762, 9.378658296173002

        //48.13692°N 11.98783°E 32 zona 
        //48.13705°N 12.00249°E 33 zona
        float UpperLat = 0.0009978417767449134f;
        float UpperLong = 9.374125157568141f;
        float LowerLat = -0.006765411712556136f;
        float LowerLong = 9.367515224307253f;

        Map map = new Map();

        map.CalculateMapWidthAndHeight(UpperLat, UpperLong, LowerLat, LowerLong);

        Debug.Log(map.Width);
        Debug.Log(map.Height);

        ObjectMap objectMap = new ObjectMap();
        objectMap.CreateObjectMap(map);

        Debug.Log(objectMap.TranslationOfLeftLowerEdgeX);
        Debug.Log(objectMap.TranslationOfLeftLowerEdgeY);

    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
