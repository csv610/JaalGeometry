/****************************************************************************
** Meta object code from reading C++ file 'MagnifyLensDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MagnifyLensDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MagnifyLensDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMagnifyingLensDialog_t {
    QByteArrayData data[8];
    char stringdata0[95];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMagnifyingLensDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMagnifyingLensDialog_t qt_meta_stringdata_JMagnifyingLensDialog = {
    {
        QT_MOC_LITERAL(0, 0, 21), // "JMagnifyingLensDialog"
        QT_MOC_LITERAL(1, 22, 12), // "updateCamera"
        QT_MOC_LITERAL(2, 35, 0), // ""
        QT_MOC_LITERAL(3, 36, 10), // "loadCamera"
        QT_MOC_LITERAL(4, 47, 10), // "addNewLens"
        QT_MOC_LITERAL(5, 58, 13), // "captureRegion"
        QT_MOC_LITERAL(6, 72, 10), // "deleteLens"
        QT_MOC_LITERAL(7, 83, 11) // "closeDialog"

    },
    "JMagnifyingLensDialog\0updateCamera\0\0"
    "loadCamera\0addNewLens\0captureRegion\0"
    "deleteLens\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMagnifyingLensDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    6,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   44,    2, 0x08 /* Private */,
    3,    0,   45,    2, 0x08 /* Private */,
    4,    0,   46,    2, 0x08 /* Private */,
    5,    0,   47,    2, 0x08 /* Private */,
    6,    0,   48,    2, 0x08 /* Private */,
    7,    0,   49,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JMagnifyingLensDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMagnifyingLensDialog *_t = static_cast<JMagnifyingLensDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->updateCamera();
            break;
        case 1:
            _t->loadCamera();
            break;
        case 2:
            _t->addNewLens();
            break;
        case 3:
            _t->captureRegion();
            break;
        case 4:
            _t->deleteLens();
            break;
        case 5:
            _t->closeDialog();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMagnifyingLensDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMagnifyingLensDialog.data,
        qt_meta_data_JMagnifyingLensDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMagnifyingLensDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMagnifyingLensDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMagnifyingLensDialog.stringdata0))
        return static_cast<void*>(const_cast< JMagnifyingLensDialog*>(this));
    if (!strcmp(_clname, "Ui::MagnifyingLensDialog"))
        return static_cast< Ui::MagnifyingLensDialog*>(const_cast< JMagnifyingLensDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMagnifyingLensDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
