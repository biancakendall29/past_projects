#include "helpers.h"
#include <math.h>

// Convert image to grayscale
void grayscale(int height, int width, RGBTRIPLE image[height][width])
{
    float average = 0;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            average = (image[i][j].rgbtRed + image[i][j].rgbtGreen + image[i][j].rgbtBlue) / 3.0; //calculate average of three colours
            image[i][j].rgbtRed = round(average); //assign averages, rounded to the nearest integer
            image[i][j].rgbtGreen = round(average);
            image[i][j].rgbtBlue = round(average);
        }
    }
    return;
}

// Convert image to sepia
void sepia(int height, int width, RGBTRIPLE image[height][width])
{
    float sepiaRed, sepiaGreen, sepiaBlue = 0;
    int originalRed, originalGreen, originalBlue = 0;

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            originalRed = image[i][j].rgbtRed;
            originalGreen = image[i][j].rgbtGreen;
            originalBlue = image[i][j].rgbtBlue;

            sepiaRed = (0.393 * originalRed) + (0.769 * originalGreen) + (0.189 * originalBlue); //formulas to convert to sepia
            sepiaGreen = (0.349 * originalRed) + (0.686 * originalGreen) + (0.168 * originalBlue);
            sepiaBlue = (0.272 * originalRed) + (0.534 * originalGreen) + (0.131 * originalBlue);

            if (sepiaRed > 255) //capping the pixel colour value at the maximum of 255
            {
                sepiaRed = 255;
            }

            if (sepiaGreen > 255)
            {
                sepiaGreen = 255;
            }

            if (sepiaBlue > 255)
            {
                sepiaBlue = 255;
            }

            image[i][j].rgbtRed = round(sepiaRed); //assigning new sepia colour to original images pixels
            image[i][j].rgbtGreen = round(sepiaGreen);
            image[i][j].rgbtBlue = round(sepiaBlue);
        }
    }
    return;
}

// Reflect image horizontally
void reflect(int height, int width, RGBTRIPLE image[height][width])
{
    RGBTRIPLE swap;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width / 2; j++)
        {
            swap = image[i][j]; //swap pixels from first column with last column and
            image[i][j] = image[i][width - 1 - j]; //second column with second last etc, and
            image[i][width - 1 - j] = swap; //iterate until half the width of the image
        }
    }
    return;
}

// Blur image
void blur(int height, int width, RGBTRIPLE image[height][width])
{
    int R1red, R1green, R1blue, R2red, R2green, R2blue, R3red, R3green, R3blue = 0;
    RGBTRIPLE copy[height][width];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            copy[i][j] = image[i][j]; //copy pixels from original image to a copy
            if (i > 0 && i < height - 1 && j > 0 && j < width - 1) //for pixels in middle of image
            {
                R1red = image[i - 1][j - 1].rgbtRed + image[i - 1][j].rgbtRed + image[i - 1][j + 1].rgbtRed;
                R2red = image[i][j - 1].rgbtRed + image[i][j].rgbtRed + image[i][j + 1].rgbtRed;
                R3red = image[i + 1][j - 1].rgbtRed + image[i + 1][j].rgbtRed + image[i + 1][j + 1].rgbtRed;

                R1green = image[i - 1][j - 1].rgbtGreen + image[i - 1][j].rgbtGreen + image[i - 1][j + 1].rgbtGreen;
                R2green = image[i][j - 1].rgbtGreen + image[i][j].rgbtGreen + image[i][j + 1].rgbtGreen;
                R3green = image[i + 1][j - 1].rgbtGreen + image[i + 1][j].rgbtGreen + image[i + 1][j + 1].rgbtGreen;

                R1blue = image[i - 1][j - 1].rgbtBlue + image[i - 1][j].rgbtBlue + image[i - 1][j + 1].rgbtBlue;
                R2blue = image[i][j - 1].rgbtBlue + image[i][j].rgbtBlue + image[i][j + 1].rgbtBlue;
                R3blue = image[i + 1][j - 1].rgbtBlue + image[i + 1][j].rgbtBlue + image[i + 1][j + 1].rgbtBlue;

                copy[i][j].rgbtRed = round((R1red + R2red + R3red) / 9.0);
                copy[i][j].rgbtGreen = round((R1green + R2green + R3green) / 9.0);
                copy[i][j].rgbtBlue = round((R1blue + R2blue + R3blue) / 9.0);
            }
            if (i == 0 && j > 0 && j < width - 1) //for pixels on top row of image
            {
                R2red = image[i][j - 1].rgbtRed + image[i][j].rgbtRed + image[i][j + 1].rgbtRed;
                R3red = image[i + 1][j - 1].rgbtRed + image[i + 1][j].rgbtRed + image[i + 1][j + 1].rgbtRed;

                R2green = image[i][j - 1].rgbtGreen + image[i][j].rgbtGreen + image[i][j + 1].rgbtGreen;
                R3green = image[i + 1][j - 1].rgbtGreen + image[i + 1][j].rgbtGreen + image[i + 1][j + 1].rgbtGreen;

                R2blue = image[i][j - 1].rgbtBlue + image[i][j].rgbtBlue + image[i][j + 1].rgbtBlue;
                R3blue = image[i + 1][j - 1].rgbtBlue + image[i + 1][j].rgbtBlue + image[i + 1][j + 1].rgbtBlue;

                copy[i][j].rgbtRed = round((R2red + R3red) / 6.0);
                copy[i][j].rgbtGreen = round((R2green + R3green) / 6.0);
                copy[i][j].rgbtBlue = round((R2blue + R3blue) / 6.0);
            }
            if (i == height - 1 && j > 0 && j < width - 1) //for pixels on bottom row of image
            {
                R1red = image[i - 1][j - 1].rgbtRed + image[i - 1][j].rgbtRed + image[i - 1][j + 1].rgbtRed;
                R2red = image[i][j - 1].rgbtRed + image[i][j].rgbtRed + image[i][j + 1].rgbtRed;

                R1green = image[i - 1][j - 1].rgbtGreen + image[i - 1][j].rgbtGreen + image[i - 1][j + 1].rgbtGreen;
                R2green = image[i][j - 1].rgbtGreen + image[i][j].rgbtGreen + image[i][j + 1].rgbtGreen;

                R1blue = image[i - 1][j - 1].rgbtBlue + image[i - 1][j].rgbtBlue + image[i - 1][j + 1].rgbtBlue;
                R2blue = image[i][j - 1].rgbtBlue + image[i][j].rgbtBlue + image[i][j + 1].rgbtBlue;

                copy[i][j].rgbtRed = round((R1red + R2red) / 6.0);
                copy[i][j].rgbtGreen = round((R1green + R2green) / 6.0);
                copy[i][j].rgbtBlue = round((R1blue + R2blue) / 6.0);
            }
            if (j == 0 && i > 0 && i < height - 1) //for pixels in first column of image
            {
                R1red = image[i - 1][j].rgbtRed + image[i - 1][j + 1].rgbtRed;
                R2red = image[i][j].rgbtRed + image[i][j + 1].rgbtRed;
                R3red = image[i + 1][j].rgbtRed + image[i + 1][j + 1].rgbtRed;

                R1green = image[i - 1][j].rgbtGreen + image[i - 1][j + 1].rgbtGreen;
                R2green = image[i][j].rgbtGreen + image[i][j + 1].rgbtGreen;
                R3green = image[i + 1][j].rgbtGreen + image[i + 1][j + 1].rgbtGreen;

                R1blue = image[i - 1][j].rgbtBlue + image[i - 1][j + 1].rgbtBlue;
                R2blue = image[i][j].rgbtBlue + image[i][j + 1].rgbtBlue;
                R3blue = image[i + 1][j].rgbtBlue + image[i + 1][j + 1].rgbtBlue;

                copy[i][j].rgbtRed = round((R1red + R2red + R3red) / 6.0);
                copy[i][j].rgbtGreen = round((R1green + R2green + R3green) / 6.0);
                copy[i][j].rgbtBlue = round((R1blue + R2blue + R3blue) / 6.0);
            }
            if (j == width - 1 && i > 0 && i < height - 1) //for pixels in last column of image
            {
                R1red = image[i - 1][j - 1].rgbtRed + image[i - 1][j].rgbtRed;
                R2red = image[i][j - 1].rgbtRed + image[i][j].rgbtRed;
                R3red = image[i + 1][j - 1].rgbtRed + image[i + 1][j].rgbtRed;

                R1green = image[i - 1][j - 1].rgbtGreen + image[i - 1][j].rgbtGreen;
                R2green = image[i][j - 1].rgbtGreen + image[i][j].rgbtGreen;
                R3green = image[i + 1][j - 1].rgbtGreen + image[i + 1][j].rgbtGreen;

                R1blue = image[i - 1][j - 1].rgbtBlue + image[i - 1][j].rgbtBlue;
                R2blue = image[i][j - 1].rgbtBlue + image[i][j].rgbtBlue;
                R3blue = image[i + 1][j - 1].rgbtBlue + image[i + 1][j].rgbtBlue;

                copy[i][j].rgbtRed = round((R1red + R2red + R3red) / 6.0);
                copy[i][j].rgbtGreen = round((R1green + R2green + R3green) / 6.0);
                copy[i][j].rgbtBlue = round((R1blue + R2blue + R3blue) / 6.0);
            }
            if (i == 0 && j == 0) //for pixels in top left corner of image
            {
                R2red = image[i][j].rgbtRed + image[i][j + 1].rgbtRed;
                R3red = image[i + 1][j].rgbtRed + image[i + 1][j + 1].rgbtRed;

                R2green = image[i][j].rgbtGreen + image[i][j + 1].rgbtGreen;
                R3green = image[i + 1][j].rgbtGreen + image[i + 1][j + 1].rgbtGreen;

                R2blue = image[i][j].rgbtBlue + image[i][j + 1].rgbtBlue;
                R3blue = image[i + 1][j].rgbtBlue + image[i + 1][j + 1].rgbtBlue;

                copy[i][j].rgbtRed = round((R2red + R3red) / 4.0);
                copy[i][j].rgbtGreen = round((R2green + R3green) / 4.0);
                copy[i][j].rgbtBlue = round((R2blue + R3blue) / 4.0);
            }
            if (i == 0 && j == width - 1) //for pixels in top right corner of image
            {
                R2red = image[i][j - 1].rgbtRed + image[i][j].rgbtRed;
                R3red = image[i + 1][j - 1].rgbtRed + image[i + 1][j].rgbtRed;

                R2green = image[i][j - 1].rgbtGreen + image[i][j].rgbtGreen;
                R3green = image[i + 1][j - 1].rgbtGreen + image[i + 1][j].rgbtGreen;

                R2blue = image[i][j - 1].rgbtBlue + image[i][j].rgbtBlue;
                R3blue = image[i + 1][j - 1].rgbtBlue + image[i + 1][j].rgbtBlue;

                copy[i][j].rgbtRed = round((R2red + R3red) / 4.0);
                copy[i][j].rgbtGreen = round((R2green + R3green) / 4.0);
                copy[i][j].rgbtBlue = round((R2blue + R3blue) / 4.0);
            }
            if (i == height - 1 && j == 0) //for pixels in bottem left corner of image
            {
                R1red = image[i - 1][j].rgbtRed + image[i - 1][j + 1].rgbtRed;
                R2red = image[i][j].rgbtRed + image[i][j + 1].rgbtRed;

                R1green = image[i - 1][j].rgbtGreen + image[i - 1][j + 1].rgbtGreen;
                R2green = image[i][j].rgbtGreen + image[i][j + 1].rgbtGreen;

                R1blue = image[i - 1][j].rgbtBlue + image[i - 1][j + 1].rgbtBlue;
                R2blue = image[i][j].rgbtBlue + image[i][j + 1].rgbtBlue;

                copy[i][j].rgbtRed = round((R1red + R2red) / 4.0);
                copy[i][j].rgbtGreen = round((R1green + R2green) / 4.0);
                copy[i][j].rgbtBlue = round((R1blue + R2blue) / 4.0);
            }
            if (i == height - 1 && j == width - 1) //for pixels in bottom right corner of image
            {
                R1red = image[i - 1][j - 1].rgbtRed + image[i - 1][j].rgbtRed;
                R2red = image[i][j - 1].rgbtRed + image[i][j].rgbtRed;

                R1green = image[i - 1][j - 1].rgbtGreen + image[i - 1][j].rgbtGreen;
                R2green = image[i][j - 1].rgbtGreen + image[i][j].rgbtGreen;

                R1blue = image[i - 1][j - 1].rgbtBlue + image[i - 1][j].rgbtBlue;
                R2blue = image[i][j - 1].rgbtBlue + image[i][j].rgbtBlue;

                copy[i][j].rgbtRed = round((R1red + R2red) / 4.0);
                copy[i][j].rgbtGreen = round((R1green + R2green) / 4.0);
                copy[i][j].rgbtBlue = round((R1blue + R2blue) / 4.0);
            }
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            image[i][j].rgbtRed = copy[i][j].rgbtRed; //copy new values to original array
            image[i][j].rgbtGreen = copy[i][j].rgbtGreen;
            image[i][j].rgbtBlue = copy[i][j].rgbtBlue;
        }
    }
    return;
}
