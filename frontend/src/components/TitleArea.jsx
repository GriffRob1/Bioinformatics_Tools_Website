import React from 'react';
import dna_code_image from '../images/dna-code-image.jpg'
import { useNavigate } from "react-router-dom";

const sectionStyle = {
    backgroundImage: 'url(' + dna_code_image + ')'
}

export default function TitleArea() {

    const navigate = useNavigate();
    const navigateToAbout = () =>{
        navigate('/about');
    }

    return (
        <section className='container title-container' style={sectionStyle}>
            <div className='title-area'>
                <h1 className='title'>Bioinformatics Tools</h1>
                <p>Essential tools for analysing biological data</p>
                <button onClick={navigateToAbout}>What is Bioinformatics?</button>
            </div>
        </section>
    )
}