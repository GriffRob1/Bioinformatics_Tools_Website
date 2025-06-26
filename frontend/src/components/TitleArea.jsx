import '../styles/TitleArea.css';
import React from 'react';
import BlueButton from "./BlueButton";

const sectionStyle = {
    backgroundImage: 'url(/images/dna-code-image.jpg)'
}

export default function TitleArea() {

    return (
        <section className={'container title-container'} style={sectionStyle}>
            <div className={'title-area'}>
                <h1 className={'title'}>Bioinformatics Tools</h1>
                <p>Essential tools for analysing biological data</p>
                <BlueButton URL={'/about'} buttonClass={'what-is-bioinformatics-button'}>What is Bioinformatics?</BlueButton>
            </div>
        </section>
    )
}